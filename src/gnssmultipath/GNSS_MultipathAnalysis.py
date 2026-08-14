"""
This is the main module for running the software GNSS_Multipath_Analysis.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""

import os
import warnings
import logging
import time
from typing import Union, List
import numpy as np
from tqdm import tqdm
from gnssmultipath.readers.readRinexObs import readRinexObs
from gnssmultipath.constants import SYSTEM_BANDS as _SYSTEM_BANDS_BY_CODE, build_frequency_overview
from gnssmultipath.Geodetic_functions import gpstime_to_utc_datefmt, gpstime2date
from gnssmultipath.computeSatElevAzimuth_fromNav import computeSatElevAzimuth_fromNav
from gnssmultipath.signalAnalysis import SignalAnalyzer
from gnssmultipath.detectClockJumps import detectClockJumps
from gnssmultipath.utils.writeOutputFile import writeOutputFile
from gnssmultipath.utils.createCSVfile import createCSVfile
from gnssmultipath.plot.make_polarplot import make_polarplot, make_skyplot, make_polarplot_SNR, plot_SNR_wrt_elev
from gnssmultipath.plot.make_polarplot import make_polarplot_dont_use_TEX, make_skyplot_dont_use_TEX, make_polarplot_SNR_dont_use_TEX, plot_SNR_wrt_elev_dont_use_TEX
from gnssmultipath.plot.plotResults import plotResults, plotResults_dont_use_TEX, make_barplot, make_barplot_dont_use_TEX
from gnssmultipath.plot.SkyPlotSummary import make_skyplot_summary
from gnssmultipath.utils.PickleHandler import PickleHandler
from gnssmultipath.PreciseSatCoords import PreciseSatCoords
from gnssmultipath.SP3PositionEstimator import SP3PositionEstimator

warnings.filterwarnings("ignore")


# ── Module-level constants ────────────────────────────────────────────────────

_CODE2NAME = {'G': 'GPS', 'R': 'GLONASS', 'E': 'Galileo', 'C': 'BeiDou'}
_NAME2CODE = {v: k for k, v in _CODE2NAME.items()}

_SYSTEM_BANDS = {
    _CODE2NAME[code]: list(bands) for code, bands in _SYSTEM_BANDS_BY_CODE.items()
}


# ── Helpers ───────────────────────────────────────────────────────────────────

def _plot_with_fallback(use_latex, latex_func, fallback_func, *args):
    """Call *latex_func*; on failure fall back to *fallback_func*. Returns False if LaTeX failed."""
    if use_latex:
        try:
            latex_func(*args)
            return True
        except Exception:
            fallback_func(*args)
            return False
    fallback_func(*args)
    return True


def compute_processing_time(start_time, end_time):
    """Computes and prints the processing time."""
    total = end_time - start_time
    h, m, s = int(total // 3600), int((total % 3600) // 60), int(total % 60)
    print(f"INFO: Finished! Processing time: {h:02d}:{m:02d}:{s:02d}")


def ismember(list_, code):
    """Return the index of *code* in *list_*, or empty list if not found."""
    indx = [idx for idx, val in enumerate(list_) if val == code]
    return indx[0] if indx else indx


def filter_common_gnss_keys(sat_pos, GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, GNSSsystems):
    """
    Filters the GNSS-related dictionaries to include only common systems present in *sat_pos* keys.

    Parameters:
    ----------
    - sat_pos (dict): Dictionary with satellite position data.
    - GNSS_obs (dict): Dictionary of GNSS observations.
    - GNSS_LLI (dict): Dictionary of GNSS Loss of Lock Indicators.
    - GNSS_SS (dict): Dictionary of GNSS Signal Strengths.
    - GNSS_SVs (dict): Dictionary of GNSS Space Vehicles.
    - GNSSsystems (dict): Dictionary of GNSS systems.

    Returns:
    --------
    - Tuple[dict, dict, dict, dict, dict]: Filtered GNSS dictionaries with common systems.
    """
    common = set(sat_pos.keys()) & set(GNSSsystems.values())
    return (
        {k: v for k, v in GNSS_obs.items() if k in common},
        {k: v for k, v in GNSS_LLI.items() if k in common},
        {k: v for k, v in GNSS_SS.items() if k in common},
        {k: v for k, v in GNSS_SVs.items() if k in common},
        {k: v for k, v in GNSSsystems.items() if v in common},
    )


# ── Main function ─────────────────────────────────────────────────────────────

def GNSS_MultipathAnalysis(rinObsFilename: str,
                          broadcastNav1: Union[str, None] = None,
                          broadcastNav2: Union[str, None] = None,
                          broadcastNav3: Union[str, None] = None,
                          broadcastNav4: Union[str, None] = None,
                          sp3NavFilename_1: Union[str, None] = None,
                          sp3NavFilename_2: Union[str, None] = None,
                          sp3NavFilename_3: Union[str, None] = None,
                          desiredGNSSsystems: Union[List[str], None] = None,
                          phaseCodeLimit: Union[float, int, None] = None,
                          ionLimit: Union[float, None] = None,
                          cutoff_elevation_angle: Union[int, None] = None,
                          outputDir: Union[str, None] = None,
                          plotEstimates: bool = True,
                          plot_polarplot: bool = True,
                          include_SNR: bool = True,
                          save_results_as_pickle: bool = True,
                          save_results_as_compressed_pickle: bool = False,
                          write_results_to_csv: bool = True,
                          output_csv_delimiter: str = ';',
                          nav_data_rate: int = 60,
                          desired_obs_data_rate: Union[float, int, None] = None,
                          includeResultSummary: Union[bool, None] = None,
                          includeCompactSummary: Union[bool, None] = None,
                          includeObservationOverview: Union[bool, None] = None,
                          includeLLIOverview: Union[bool, None] = None,
                          use_LaTex: bool = True
                          ) -> dict:

    """
    GNSS Multipath Analysis
    ------------------------
    Made by: Per Helge Aarnes
    E-mail: per.helge.aarnes@gmail.com

    GNSS_MultipathAnalysis is a software for analyzing the multipath effect on Global Navigation Satellite Systems (GNSS) and
    is based on the MATLAB software "GNSS_Receiver_QC_2020" made by Bjørn-Eirik Roald.
    This is the main function of the software that, through the help of other functions:

      - reads RINEX observation files
      - reads SP3 satellite navigation files (if a SP3 file is fed in)
      - reads rinex navigation files (if a rinex navigation file is fed in)
      - computes elevation angles of satellites for every epoch
      - makes estimates of multipath, ionospheric delay, and cycle slips
          for all signals in RINEX observation file
      - plots estimates if user choose to do it
      - computes and stores statistics on estimates
      - writes an output files containing results

    This function calls on a range of functions. These in turn call on
    further more functions. The function called upon directly in
    GNSS_MultipathAnalysis are:

      - readRinexObs.py
      - computeSatElevations.py
      - computeSatElevAzimuth_fromNav.py
      - signalAnalysis.py
      - plotEstimates.py
      - make_barplot.py
      - make_polarplot.py
      - writeOutputFile.py

    --------------------------------------------------------------------------------------------------------------------------

    INPUTS:
    ------

    rinObsFilename:           string. Path to RINEX 3 observation file

    sp3NavFilename_1:         string. Path to first SP3 navigation file

    sp3NavFilename_2:         string. Path to second SP3 navigation file (optional).

    sp3NavFilename_3:         string. Path to third SP3 navigation file (optional).

    desiredGNSSsystems:       List with the desired GNSS system. Ex ['G','R'] if you want to
                              only run the analysis on GPS and GLONASS. Default: All systems. (if set to None) (optional)

    phaseCodeLimit:           critical limit that indicates cycle slip for
                              phase-code combination. Unit: m/s. If set to 0,
                              default value of 6.667 m/s will be used (optional)

    ionLimit:                 critical limit that indicates cycle slip for
                              the rate of change of the ionopheric delay.
                              Unit: m/s. If set to 0, default value of 0.0667 m/s will be used (optional)

    cutoff_elevation_angle    Critical cutoff angle for satellite elevation angles, degrees
                              Estimates where satellite elevation angle
                              is lower than cutoff are removed, so are
                              estimated slip periods (optional)

    outputDir:                string. Path to directory where output file
                              should be generated. If user does not wish to
                              specify output directory, this variable should
                              be empty string, "". In this case the output file
                              will be generated in sub-directory inside same
                              directory as GNSS_Receiver_QC_2020.m (optional)

    plotEstimates:            boolean. False if user desires estimates not to be
                              ploted. True by default. (optional)

    plot_polarplot:           boolean. True or False. If not defined polarplots will be made (optional)


    include_SNR:              boolean. If not defined, SNR from Rinex obs file will NOT be used (optional)

    save_results_as_pickle:   boolean. If True, the results will be stored as dictionary in form of a binary pickle file. Default set to True.


    save_results_as_compressed_pickle : boolean. If True, the results will be stored as dictionary in form of a binary compressed pickle file (zstd compression). Default set to False.

    write_results_to_csv: boolean. If True, a subset of the results will be exported as a CSV file. Default is True.

    output_csv_delimiter:     str. Set the delimiter of the CSV file. Default is semi colon (;).


    nav_data_rate:            integer. The desired data rate of ephemerides given in minutes. Default is 60 min. The purpose with this
                              parameter is to speed up processing time. Both reading the RINEX navigation file and looping through the
                              ephemerides aftwerward will be significatnly faster by increasing this value. Note: A too high value will
                              degrade the accuracy of the interploated satellite coordinates.

    desired_obs_data_rate:    float | int | None. Desired observation data rate in **seconds** for down-sampling (decimating)
                              the RINEX observation file at read-time. Use this to reduce the amount of data when the file
                              has a higher native rate than required (e.g. read a 1 Hz file at 30 s by passing ``30``).
                              The reader keeps every ``round(desired_obs_data_rate / native_tInterval)``-th epoch starting
                              from the first epoch. ``None`` (default) keeps all epochs. A value smaller than or equal to the
                              file's native interval also keeps all epochs (no up-sampling is performed).

    includeResultSummary:     boolean. 1 if user desires output file to
                              include more detailed overview of statistics,
                              including for each individual satellites.
                              0 otherwise (optional)

    includeCompactSummary:    boolean. 1 if user desired output file to
                              include more compact overview og statistics. (optional)

    includeObservationOverview: boolean. 1 if user desires output file to
                                include overview of obseration types observed
                                by each satellite. 0 otherwise (optional)

    use_LaTex:                 boolean. Will use TeX as an interpreter in plots. Default set to true. "Requires TeX installed on computer".

    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:

    analysisResults:          A dictionary that contains alls results of all analysises, for all GNSS systems.


    The software is also returning results file. A report provided as a text file, and a CSV file with the estimated values.
    --------------------------------------------------------------------------------------------------------------------------
    """
    start_time = time.time()

    # ── Input validation ──────────────────────────────────────────────────────
    if broadcastNav1 is None and sp3NavFilename_1 is None:
        raise RuntimeError("No SP3 or navigation file is defined! This is "
                           "mandatory for this software, so please add one of them.")

    if broadcastNav1 is not None and sp3NavFilename_1 is not None:
        raise RuntimeError("You defined both a navigation file and a SP3 file. Please "
                           "choose between using broadcast ephemerides or precise.")

    for name, val, allow_none in [
        ('rinObsFilename',  rinObsFilename,  False),
        ('sp3NavFilename_1', sp3NavFilename_1, True),
        ('sp3NavFilename_2', sp3NavFilename_2, True),
        ('sp3NavFilename_3', sp3NavFilename_3, True),
    ]:
        if allow_none and val is None:
            continue
        if not isinstance(val, str):
            print(f'ERROR(GNSS_MultipathAnalysis): The input variable {name} must be a string\n'
                  f'Argument is now of type {type(val)}\n')
            return

    if not os.path.isfile(rinObsFilename):
        print('ERROR(GNSS_MultipathAnalysis): RINEX observation file can not be found. '
              'Please check that the path is correct.\n')
        return

    # ── Apply defaults ────────────────────────────────────────────────────────
    broadcastNav1      = broadcastNav1 or ""
    broadcastNav2      = broadcastNav2 or ""
    broadcastNav3      = broadcastNav3 or ""
    broadcastNav4      = broadcastNav4 or ""
    sp3NavFilename_1   = sp3NavFilename_1 or ""
    sp3NavFilename_2   = sp3NavFilename_2 or ""
    sp3NavFilename_3   = sp3NavFilename_3 or ""
    phaseCodeLimit     = phaseCodeLimit if phaseCodeLimit is not None else 4 / 60 * 100
    ionLimit           = ionLimit if ionLimit is not None else 4 / 60
    cutoff_elevation_angle = cutoff_elevation_angle if cutoff_elevation_angle is not None else 0
    outputDir          = outputDir or 'Output_Files'
    includeResultSummary       = includeResultSummary if includeResultSummary is not None else 1
    includeCompactSummary      = includeCompactSummary if includeCompactSummary is not None else 1
    includeObservationOverview = includeObservationOverview if includeObservationOverview is not None else 1
    includeLLIOverview         = includeLLIOverview if includeLLIOverview is not None else 1

    includeAllGNSSsystems = 1 if desiredGNSSsystems is None else 0
    desiredGNSSsystems    = desiredGNSSsystems or ["G", "R", "E", "C"]

    if save_results_as_pickle and save_results_as_compressed_pickle:
        save_results_as_pickle = False

    use_sp3 = bool(sp3NavFilename_1)

    # ── Warn about missing optional SP3 files ─────────────────────────────────
    for label, path in [("First", sp3NavFilename_1), ("Second", sp3NavFilename_2), ("Third", sp3NavFilename_3)]:
        if path and not os.path.isfile(path):
            print(f'WARNING: {label} SP3 Navigation file can not be found.\n')

    # ── Setup directories & logging ───────────────────────────────────────────
    os.makedirs(outputDir, exist_ok=True)
    graphDir = os.path.join(outputDir, 'Graphs')
    os.makedirs(graphDir, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(outputDir, 'Logfile.log'),
        filemode='w', level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    # ── Read RINEX observation file ───────────────────────────────────────────
    desiredObsCodes = ["C", "L", "S"] if include_SNR else ["C", "L"]

    rinex_data = readRinexObs(rinObsFilename, readSS=1, readLLI=1,
                     includeAllGNSSsystems=includeAllGNSSsystems,
                     includeAllObsCodes=0,
                     desiredGNSSsystems=desiredGNSSsystems,
                     desiredObsCodes=desiredObsCodes,
                     desiredObsBands=list(range(1, 10)),
                     desired_data_rate=desired_obs_data_rate)
    GNSS_obs = rinex_data.GNSS_obs
    GNSS_LLI = rinex_data.GNSS_LLI
    GNSS_SS = rinex_data.GNSS_SS
    GNSS_SVs = rinex_data.GNSS_SVs
    time_epochs = rinex_data.time_epochs
    nepochs = rinex_data.nepochs
    GNSSsystems = rinex_data.GNSSsystems
    obsCodes = rinex_data.obsCodes
    approxPosition = rinex_data.approxPosition
    max_sat = rinex_data.max_sat
    tInterval = rinex_data.tInterval
    markerName = rinex_data.markerName
    rinexVersion = rinex_data.rinexVersion
    recType = rinex_data.recType
    timeSystem = rinex_data.timeSystem
    leapSec = rinex_data.leapSec
    gnssType = rinex_data.gnssType
    rinexProgr = rinex_data.rinexProgr
    rinexDate = rinex_data.rinexDate
    antDelta = rinex_data.antDelta
    tFirstObs = rinex_data.tFirstObs
    tLastObs = rinex_data.tLastObs
    clockOffsetsON = rinex_data.clockOffsetsON
    GLO_Slot2ChannelMap = rinex_data.GLO_Slot2ChannelMap
    success = rinex_data.success

    # ── Compute satellite positions & elevation angles ────────────────────────
    sat_pos = {}
    estimated_position, stats = None, None
    latex_installed = True
    glo_fcn = None

    if use_sp3:
        x_rec_approx, y_rec_approx, z_rec_approx = np.atleast_2d(approxPosition).flatten()
        sp3_files = [f for f in [sp3NavFilename_1, sp3NavFilename_2, sp3NavFilename_3] if f]

        sat_obj = PreciseSatCoords(sp3_files, time_epochs=time_epochs, GNSSsystems=GNSSsystems)
        df_sat_coordinates = sat_obj.satcoords

        if all(c == 0 for c in [x_rec_approx, y_rec_approx, z_rec_approx]):
            desired_time = np.array(gpstime2date(time_epochs[0, 0], round(time_epochs[0, 1], 6)))
            position_estimator = SP3PositionEstimator(
                df_sat_coordinates, desired_time=desired_time,
                GNSS_obs=GNSS_obs, time_epochs=time_epochs,
                GNSSsystems=GNSSsystems, obsCodes=obsCodes,
                sp3_metadata_dict=sat_obj.sp3_metadata_dict)
            estimated_position, stats = position_estimator.estimate_position()
            x_rec_approx, y_rec_approx, z_rec_approx, _ = estimated_position.flatten()

        df_az_el = sat_obj.compute_azimuth_and_elevation(
            receiver_position=(x_rec_approx, y_rec_approx, z_rec_approx), drop_below_horizon=True)
        sat_dict = sat_obj.create_satellite_data_dict(df_sat_coordinates, df_az_el)

        sat_coordinates, sat_elevation_angles, sat_azimut_angles = {}, {}, {}
        for idx, (system, data) in enumerate(sat_dict.items()):
            sat_coordinates[system] = data.get('coordinates', {})
            sat_elevation_angles[idx] = data.get('elevation', None)
            sat_azimut_angles[idx] = data.get('azimuth', None)

        # Free intermediate per-satellite DataFrames now that the dense
        # NumPy arrays have been extracted.
        del sat_dict, df_az_el, df_sat_coordinates

    else:
        nav_files = [broadcastNav1, broadcastNav2, broadcastNav3, broadcastNav4]
        sat_pos, glo_fcn, estimated_position, stats = computeSatElevAzimuth_fromNav(
            nav_files, approxPosition, GNSS_SVs, time_epochs, nav_data_rate,
            GNSS_obs, GNSSsystems, obsCodes)

        sat_elevation_angles = {}
        sat_pos_dummy = sat_pos.copy()

        GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, GNSSsystems = filter_common_gnss_keys(
            sat_pos_dummy, GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, GNSSsystems)

        for sys in range(len(GNSSsystems)):
            sys_code = GNSSsystems[sys + 1]
            if sys_code not in sat_pos_dummy:
                sys_name = _CODE2NAME[sys_code]
                available = [_CODE2NAME[s] for s in sat_pos_dummy]
                raise KeyError(f'GNSS system "{sys_name}" is not present in the RINEX navigation file: {available}')

            sat_elevation_angles[sys] = (sat_pos_dummy[sys_code]['elevation'][:, 0:37]
                                         if sys_code != 'C'
                                         else sat_pos_dummy[sys_code]['elevation'])

        missing_sys = []
        for key, sys_code in list(GNSSsystems.items()):
            if len(sat_pos[sys_code]['position']) == 0:
                del sat_pos[sys_code], GNSSsystems[key]
                missing_sys.append(sys_code)
        if missing_sys:
            print(f'\n\nSystems {missing_sys} does not exist in navigation file! '
                  '\nMultipath analysis for these systems is therefore not possible. '
                  '\nConsider using another navigation file.\n\n')

    # ── Build carrier frequency overview ──────────────────────────────────────
    nGNSSsystems = len(GNSSsystems)

    if "R" in GNSSsystems.values() and not isinstance(GLO_Slot2ChannelMap, dict):
        # Fall back to the k-numbers decoded from the navigation file
        GLO_Slot2ChannelMap = glo_fcn

    frequencyOverview = build_frequency_overview(GNSSsystems, GLO_Slot2ChannelMap)

    # ── Build observation code overview (per system, per band) ────────────────
    obsCodeOverview = {}
    for sys_key, sys_code in GNSSsystems.items():
        obsCodeOverview[sys_key] = {str(b): [] for b in range(1, 10)}
        codes = [c for c in obsCodes[sys_key][sys_code] if c[0] in ('C', 'P')]
        for band in {c[1] for c in codes}:
            obsCodeOverview[sys_key][band] = [c for c in codes if c[1] == band]

    # ── Initialize results dictionary ─────────────────────────────────────────
    nCodes_Total = 0
    analysisResults = {
        'nGNSSsystem': nGNSSsystems,
        'GNSSsystems': list(GNSSsystems.values()),
    }

    for sys in range(nGNSSsystems):
        sys_name = _CODE2NAME[GNSSsystems[sys + 1]]
        possible_bands = _SYSTEM_BANDS[sys_name]
        band_names = [f"Band_{b}" for b in possible_bands]

        # Observation overview per satellite
        obs_overview = {}
        for sat_num in range(1, int(max_sat[sys][0]) + 1):
            sat_entry = {'n_possible_bands': len(possible_bands), 'Bands': list(band_names)}
            for bn in band_names:
                sat_entry[bn] = ""
            obs_overview[f'Sat_{sat_num}'] = sat_entry

        # Band structure
        current_sys_dict = {'observationOverview': obs_overview}
        bands_list = []
        for band_num in range(1, 10):
            codes = obsCodeOverview[sys + 1][str(band_num)]
            if codes:
                band_name = f"Band_{band_num}"
                bands_list.append(band_name)
                nCodes_Total += len(codes)
                current_sys_dict[band_name] = {'nCodes': len(codes), 'Codes': codes}
        current_sys_dict['nBands'] = len(bands_list)
        current_sys_dict['Bands'] = bands_list
        analysisResults[sys_name] = current_sys_dict

    # ── Execute signal analysis ───────────────────────────────────────────────
    codeNum = 0
    bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'

    for sys in range(nGNSSsystems):
        currentGNSSsystem = GNSSsystems[sys + 1]
        GNSSsystemName = _CODE2NAME[currentGNSSsystem]
        current_sys_dict = analysisResults[GNSSsystemName]
        nBands = current_sys_dict['nBands']

        n_signals = sum(current_sys_dict[current_sys_dict['Bands'][b]]['nCodes'] for b in range(nBands))
        pbar = tqdm(total=n_signals,
                    desc=f'Currently processing all available signals for {GNSSsystemName}',
                    position=0, leave=True, bar_format=bar_format)

        # Precompute per-obscode "has any non-zero data" bitmap for this system.
        # Replaces a large np.stack(GNSS_obs[...].values()) (hundreds of MB at 1 Hz daily).
        # Each per-epoch array is (nsats, nobscodes); we OR-reduce across epochs into
        # a single (nobscodes,) bool vector with negligible transient memory.
        sys_obs_codes = obsCodes[sys + 1][currentGNSSsystem]
        has_data = np.zeros(len(sys_obs_codes), dtype=bool)
        for _epoch_arr in GNSS_obs[currentGNSSsystem].values():
            if has_data.all():
                break
            has_data |= np.any(_epoch_arr != 0, axis=0)

        for bandNumInd in range(nBands):
            current_band_dict = current_sys_dict[current_sys_dict['Bands'][bandNumInd]]
            currentBandName = current_sys_dict['Bands'][bandNumInd]

            for i in range(len(current_band_dict['Codes'])):
                range1_Code = current_band_dict['Codes'][i]
                if not isinstance(range1_Code, str):
                    continue
                phase1_Code = "L" + range1_Code[1:]
                codeNum += 1
                obs_codes_list = obsCodes[sys + 1][currentGNSSsystem]

                # Skip if matching phase observation not available
                if phase1_Code not in obs_codes_list:
                    pbar.update(1)
                    logger.warning(
                        f"INFO(GNSS_MultipathAnalysis): {range1_Code} code exists in RINEX observation file, "
                        f"but not {phase1_Code}\nLinear combination using this signal is not used.")
                    current_band_dict['Codes'][ismember(current_band_dict['Codes'], range1_Code)] = []
                    current_band_dict['nCodes'] -= 1
                    continue

                # Pre-compute values constant across all secondary band combinations
                range1_Code_idx = obs_codes_list.index(range1_Code)
                phase1_Code_idx = obs_codes_list.index(phase1_Code)
                current_max_sat = int(max_sat[sys].item()) if hasattr(max_sat[sys], "item") else int(max_sat[sys])

                best_nEstimates = 0
                best_currentStats = np.nan

                # Iterate through codes in other bands to find best linear combination
                for secondBandnum in range(nBands):
                    if secondBandnum == bandNumInd:
                        continue

                    other_band_dict = current_sys_dict[current_sys_dict['Bands'][secondBandnum]]

                    for k in range(len(other_band_dict['Codes'])):
                        range2_Code = other_band_dict['Codes'][k]
                        if not isinstance(range2_Code, str):
                            continue
                        phase2_Code = "L" + range2_Code[1:]

                        # Check if phase2 observation was read from RINEX observation file
                        if phase2_Code not in obs_codes_list:
                            logger.warning(
                                f"INFO(GNSS_MultipathAnalysis): {range2_Code} code exists in RINEX observation file, "
                                f"but not {phase2_Code}. Linear combinations using this signal are not used.")
                            continue

                        # Check that signals contain actual data
                        range2_Code_idx = obs_codes_list.index(range2_Code)
                        phase2_Code_idx = obs_codes_list.index(phase2_Code)
                        if not (has_data[range1_Code_idx] and has_data[phase1_Code_idx]
                                and has_data[range2_Code_idx] and has_data[phase2_Code_idx]):
                            logger.warning(
                                f"INFO(GNSS_MultipathAnalysis): One or more of the following observation codes "
                                f"{range1_Code},{phase1_Code} and {phase2_Code} ({GNSSsystemName}),"
                                " lack data for the entire observation period. "
                                "Therefore, this linear combination cannot be utilized.")
                            continue

                        # Execute the analysis of current combination of observations.
                        # NOTE: name kept distinct from the outer ``stats`` (position
                        # estimator output) so it isn't accidentally shadowed and later
                        # reused as ``Estimated_Receiver_Approx_Pos_stats``.
                        signal_stats, success = SignalAnalyzer(
                            gnss_system=currentGNSSsystem,
                            range1_code=range1_Code,
                            range2_code=range2_Code,
                            gnss_systems=GNSSsystems,
                            frequency_overview=frequencyOverview,
                            nepochs=nepochs,
                            t_interval=tInterval,
                            max_sat=current_max_sat,
                            gnss_svs=GNSS_SVs[currentGNSSsystem],
                            obs_codes=obsCodes[sys + 1],
                            gnss_obs=GNSS_obs[currentGNSSsystem],
                            gnss_lli=GNSS_LLI[currentGNSSsystem],
                            sat_elevation_angles=sat_elevation_angles[sys],
                            phase_code_limit=phaseCodeLimit,
                            ion_limit=ionLimit,
                            cutoff_elevation_angle=cutoff_elevation_angle,
                        ).run()

                        if not success:
                            return success

                        currentStats = signal_stats.to_dict()
                        current_nEstimates = currentStats['nEstimates']
                        if current_nEstimates == 0:
                            logger.warning(
                                f'INFO(GNSS_MultipathAnalysis): Estimates for signal combination '
                                f'{range1_Code}-{phase1_Code}-{phase2_Code} were not possible. '
                                'The reason could be a lack of simultaneous observations from the three signals.')

                        if current_nEstimates > best_nEstimates:
                            best_nEstimates = current_nEstimates
                            best_currentStats = currentStats

                # Store best analysis result
                current_code_dict = best_currentStats
                if not isinstance(current_code_dict, dict):
                    pbar.update(1)
                    continue

                # Update observation overview per satellite
                nSat = len(current_code_dict['range1_slip_distribution_per_sat'])
                for sat in range(nSat):
                    if current_code_dict['n_range1_obs_per_sat'][0, sat + 1] > 0:
                        sat_key = f'Sat_{sat + 1}'
                        existing = current_sys_dict['observationOverview'][sat_key][currentBandName]
                        r1_code = current_code_dict['range1_Code']
                        if existing != r1_code:
                            current_sys_dict['observationOverview'][sat_key][currentBandName] = (
                                r1_code if not existing else f"{existing}, {r1_code}")

                # Plot estimates
                if plotEstimates:
                    plot_args = (
                        current_code_dict['ion_delay_phase1'], current_code_dict['multipath_range1'],
                        current_code_dict['sat_elevation_angles'], tInterval, currentGNSSsystem,
                        current_code_dict['range1_Code'], current_code_dict['range2_Code'],
                        current_code_dict['phase1_Code'], current_code_dict['phase2_Code'], graphDir)
                    ok = _plot_with_fallback(use_LaTex, plotResults, plotResults_dont_use_TEX, *plot_args)
                    if not ok:
                        latex_installed = False
                        # Disable LaTeX for the remainder of the run to avoid
                        # retrying (and leaking memory) on every signal pair.
                        use_LaTex = False

                current_band_dict[range1_Code] = current_code_dict
                pbar.update(1)

            # Store updated band dict back into system dict
            current_sys_dict[current_sys_dict['Bands'][bandNumInd]] = current_band_dict

        # Store satellite positions for current system (SP3 path)
        if use_sp3:
            try:
                sat_pos[currentGNSSsystem] = {
                    'position':  sat_coordinates[currentGNSSsystem],
                    'azimuth':   sat_azimut_angles[sys],
                    'elevation': sat_elevation_angles[sys],
                }
            except Exception:
                pass

        analysisResults[GNSSsystemName] = current_sys_dict
        pbar.close()

    # ── Store metadata in results (loop-invariant — computed once) ────────────
    analysisResults['Sat_position'] = sat_pos

    nClockJumps, meanClockJumpInterval, stdClockJumpInterval = detectClockJumps(
        GNSS_obs, nGNSSsystems, obsCodes, time_epochs, tInterval, GNSSsystems)

    analysisResults['ExtraOutputInfo'] = {
        'rinex_obs_filename':       os.path.basename(rinObsFilename),
        'markerName':               markerName,
        'rinexVersion':             rinexVersion,
        'rinexProgr':               rinexProgr,
        'recType':                  recType,
        'Rinex_Receiver_Approx_Pos': np.atleast_2d(approxPosition).flatten().tolist(),
        'tFirstObs':                tFirstObs,
        'tLastObs':                 tLastObs,
        'tInterval':                tInterval,
        'time_epochs_gps_time':     time_epochs,
        'time_epochs_utc_time':     gpstime_to_utc_datefmt(time_epochs),
        'GLO_Slot2ChannelMap':      GLO_Slot2ChannelMap,
        'nEpochs':                  nepochs,
        'elevation_cutoff':         cutoff_elevation_angle,
        'phaseCodeLimit':           4 / 60 * 100 if phaseCodeLimit == 0 else phaseCodeLimit,
        'ionLimit':                 4 / 60 if ionLimit == 0 else ionLimit,
        'nClockJumps':              nClockJumps,
        'meanClockJumpInterval':    meanClockJumpInterval,
        'stdClockJumpInterval':     stdClockJumpInterval,
    }

    if estimated_position is not None:
        analysisResults['ExtraOutputInfo']['Estimated_Receiver_Approx_Pos'] = np.round(estimated_position.flatten(), 4).tolist()[:-1]
        analysisResults['ExtraOutputInfo']['Estimated_Receiver_Approx_Pos_stats'] = stats

    if sp3NavFilename_1:
        sp3_list = [sp3NavFilename_1, sp3NavFilename_2, sp3NavFilename_3]
        analysisResults['ExtraOutputInfo']['SP3_filename'] = [os.path.basename(f) for f in sp3_list if f]

    if broadcastNav1:
        nav_list = [broadcastNav1, broadcastNav2, broadcastNav3, broadcastNav4]
        analysisResults['ExtraOutputInfo']['rinex_nav_filename'] = [os.path.basename(f) for f in nav_list if f]

    # ── Write output file ─────────────────────────────────────────────────────
    print('\n\nINFO: Analysis complete!\n')
    if not latex_installed:
        logger.warning("INFO(GNSS_MultipathAnalysis): Use of TEX was enabled, but not installed on your computer! "
                       "Install that to get prettier text formatting in plots.")

    outputFilename = os.path.basename(rinObsFilename).split('.')[0] + '_Report.txt'
    writeOutputFile(outputFilename, outputDir, analysisResults,
                    includeResultSummary, includeCompactSummary, includeObservationOverview, includeLLIOverview)
    print(f'INFO: The output file {outputFilename} has been written.\n')

    # ── Generate plots ────────────────────────────────────────────────────────
    if plotEstimates:
        print('INFO: Making bar plot. Please wait...\n')
        _plot_with_fallback(use_LaTex, make_barplot, make_barplot_dont_use_TEX, analysisResults, graphDir)

        for sys_name in analysisResults['GNSSsystems']:
            sys_code = _NAME2CODE[sys_name]
            try:
                az = analysisResults['Sat_position'][sys_code]['azimuth']
                el = analysisResults['Sat_position'][sys_code]['elevation']
                print('INFO: Making a regular polar plot for showing azimut and elevation angle for each satellite. Please wait...')
                _plot_with_fallback(use_LaTex, make_skyplot, make_skyplot_dont_use_TEX, az, el, sys_name, graphDir)
            except Exception:
                print(f'Skyplot is not possible for {sys_name}! Missing data.')

        if plot_polarplot:
            print('INFO: Making a polar plot of the multipath effect. Please wait ...')
            _plot_with_fallback(use_LaTex, make_polarplot, make_polarplot_dont_use_TEX, analysisResults, graphDir)

    # ── SNR processing ────────────────────────────────────────────────────────
    if include_SNR:
        SNR_codes = []
        for sys_code in GNSS_obs:
            sys_key = next(k for k, v in GNSSsystems.items() if v == sys_code)
            SNR_codes = [c for c in obsCodes[sys_key][sys_code] if c.startswith('S')]

        if plotEstimates and SNR_codes:
            print('INFO: Making a plot of the Signal To Noise Ration (SNR). Please wait ...')
            snr_args = (analysisResults, GNSS_obs, GNSSsystems, obsCodes, graphDir)
            if use_LaTex:
                try:
                    make_polarplot_SNR(*snr_args)
                    plot_SNR_wrt_elev(*snr_args, tInterval)
                except Exception:
                    make_polarplot_SNR_dont_use_TEX(*snr_args)
                    plot_SNR_wrt_elev_dont_use_TEX(*snr_args, tInterval)
            else:
                make_polarplot_SNR_dont_use_TEX(*snr_args)
                plot_SNR_wrt_elev_dont_use_TEX(*snr_args, tInterval)
        else:
            logger.warning("INFO(GNSS_MultipathAnalysis): There is no SNR codes available in the RINEX files. "
                           "Plot of the Signal To Noise Ration is not possible.")

        for syst_key in GNSSsystems:
            syst_name = _CODE2NAME[GNSSsystems[syst_key]]
            analysisResults[syst_name]['SNR'] = {}
            curr_obscodes = obsCodes[syst_key][GNSSsystems[syst_key]]
            sys_obs_dict = GNSS_obs[GNSSsystems[syst_key]]
            for n, code in enumerate(curr_obscodes):
                if code.startswith('S'):
                    # Extract just the n-th obscode column from each per-epoch
                    # array and stack into (nepochs, nsats); avoids materializing
                    # the full (nepochs, nsats, nobscodes) 3D copy.
                    curr_signal = np.stack([arr[:, n] for arr in sys_obs_dict.values()])
                    analysisResults[syst_name]['SNR'][code] = np.squeeze(curr_signal)

    # ── Summary heatmaps (issue #55) ───────────────────────────────────────
    if plotEstimates:
        print('INFO: Making azimuth vs elevation summary heatmaps. Please wait...')
        try:
            make_skyplot_summary(analysisResults, graphDir, use_latex=use_LaTex)
        except Exception:
            try:
                make_skyplot_summary(analysisResults, graphDir, use_latex=False)
            except Exception as e:
                logger.warning("INFO(GNSS_MultipathAnalysis): Summary heatmaps could not be generated: %s", e)

    # ── Save results ──────────────────────────────────────────────────────────
    if write_results_to_csv:
        result_dir = os.path.join(outputDir, "Result_files_CSV")
        os.makedirs(result_dir, exist_ok=True)
        createCSV = createCSVfile(analysisResults, result_dir, output_csv_delimiter)
        createCSV.write_results_to_csv()

    if save_results_as_pickle:
        pickle_filename = 'analysisResults.pkl'
        print(f'\nINFO: The analysis results are being written to the file {pickle_filename}. Please wait..')
        PickleHandler.write_pickle(analysisResults, os.path.join(outputDir, pickle_filename))
        print(f'INFO: The analysis results has been written to the file {pickle_filename}.\n')
    elif save_results_as_compressed_pickle:
        pickle_filename = 'analysisResults.pkl.zst'
        print(f'\nINFO: The analysis results are being written to the file {pickle_filename}. Please wait..')
        PickleHandler.write_zstd_pickle(analysisResults, os.path.join(outputDir, pickle_filename))
        print(f'INFO: The analysis results has been written to the file {pickle_filename}.\n')

    end_time = time.time()
    compute_processing_time(start_time, end_time)
    logging.shutdown()

    return analysisResults


if __name__ == "__main__":
    pass
