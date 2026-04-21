"""
Module for reading RINEX observation files in v2, v3, and v4.

RINEX v4 observation files use the same record format as v3.xx,
so the v3.04 reader handles them directly. The major v4 changes
affect navigation files (new message-type headers), not observations.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com
"""
import time
import os
import re
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, Tuple

import numpy as np
from tqdm import tqdm
from gnssmultipath.Geodetic_functions import date2gpstime
from gnssmultipath.readers.GNSSObservationData import GNSSObservationData


global tFirstObs


# ── RINEX observation data container ─────────────────────────────────────────

@dataclass
class RinexObsData:
    """Container for parsed RINEX observation data (all versions).

    Consolidates the 25 return values of the legacy ``readRinexObs*``
    functions into a single, attribute-accessible object.  Backward
    compatibility is preserved via ``as_tuple()`` which returns the
    original 25-element tuple in the same order.
    """

    # Observation data
    GNSS_obs:  Dict[str, Any] = field(default_factory=dict)
    GNSS_LLI:  Any = np.nan
    GNSS_SS:   Any = np.nan
    GNSS_SVs:  Any = np.nan
    time_epochs: Any = np.nan
    nepochs: int = 0
    GNSSsystems: Dict = field(default_factory=dict)
    obsCodes: Dict = field(default_factory=dict)

    # Receiver / site metadata
    approxPosition: Any = np.nan
    max_sat: Any = np.nan
    tInterval: float = np.nan
    markerName: str = ''

    # RINEX header metadata
    rinexVersion: str = ''
    recType: str = ''
    timeSystem: str = ''
    leapSec: float = np.nan
    gnssType: str = ''
    rinexProgr: str = ''
    rinexDate: str = ''
    antDelta: Any = np.nan

    # Time bounds
    tFirstObs: Any = np.nan
    tLastObs: Any = np.nan

    # Flags and maps
    clockOffsetsON: int = 0
    GLO_Slot2ChannelMap: Any = np.nan
    success: int = 1

    # ── Tuple interop ────────────────────────────────────────────────────

    _FIELDS_ORDERED = (
        'GNSS_obs', 'GNSS_LLI', 'GNSS_SS', 'GNSS_SVs', 'time_epochs',
        'nepochs', 'GNSSsystems', 'obsCodes', 'approxPosition', 'max_sat',
        'tInterval', 'markerName', 'rinexVersion', 'recType', 'timeSystem',
        'leapSec', 'gnssType', 'rinexProgr', 'rinexDate', 'antDelta',
        'tFirstObs', 'tLastObs', 'clockOffsetsON', 'GLO_Slot2ChannelMap',
        'success',
    )

    def as_tuple(self) -> Tuple:
        """Return the legacy 25-element tuple (same order as old API)."""
        return tuple(getattr(self, f) for f in self._FIELDS_ORDERED)

    def __iter__(self):
        """Allow ``a, b, ... = readRinexObs(...)`` unpacking."""
        return iter(self.as_tuple())

    def __len__(self):
        return len(self._FIELDS_ORDERED)

    def __getitem__(self, idx):
        return self.as_tuple()[idx]

    @classmethod
    def from_tuple(cls, values: Tuple) -> 'RinexObsData':
        """Construct from an existing 25-element tuple."""
        return cls(**dict(zip(cls._FIELDS_ORDERED, values)))

    @property
    def observations(self) -> GNSSObservationData:
        """Pythonic observation accessor.

        Returns an :class:`GNSSObservationData` that provides easy,
        code-based access to the observations, LLI, and signal strength::

            obs = rinex_data.observations
            gps_c1c = obs.gps['C1C']        # ndarray [epochs, sats]
            obs.galileo.phase_codes          # ['L1X', 'L5X', ...]
        """
        if not hasattr(self, '_observation_store'):
            self._observation_store = GNSSObservationData(
                self.GNSS_obs,
                self.obsCodes,
                self.GNSSsystems,
                self.GNSS_LLI,
                self.GNSS_SS,
            )
        return self._observation_store


# ── Public dispatcher ────────────────────────────────────────────────────────

def readRinexObs(filename, readSS=None, readLLI=None,
                 includeAllGNSSsystems=None, includeAllObsCodes=None,
                 desiredGNSSsystems=None, desiredObsCodes=None,
                 desiredObsBands=None, desired_data_rate=None):
    """Read a RINEX observation file (v2 / v3 / v4).

    Dispatches to the correct version-specific reader and returns
    a :class:`RinexObsData` instance.  The object supports tuple
    unpacking, so existing call-sites that do::

        GNSS_obs, GNSS_LLI, ..., success = readRinexObs(file)

    continue to work unchanged.

    Parameters
    ----------
    desired_data_rate : float | int | None, optional
        Desired observation interval in **seconds**.  When provided, the
        epochs are decimated (down-sampled) so that the resulting interval
        is approximately ``desired_data_rate`` seconds.  Useful for
        reducing the amount of data when the file has a higher native
        rate than required (e.g. read 1 Hz data as 30 s).

        Behaviour:

        * ``None`` (default) -- keep all epochs (no decimation).
        * Value ``<=`` the file's native ``tInterval`` -- no decimation.
        * Value ``>`` ``tInterval`` -- keep every ``round(desired_data_rate /
          tInterval)``-th epoch starting from the first epoch.

        After decimation, ``RinexObsData.tInterval`` is updated to the
        effective interval and ``nepochs`` reflects the new epoch count.
    """
    if os.stat(filename).st_size == 0:
        raise ValueError('ERROR: This file seems to be empty')

    with open(filename, 'r') as fid:
        line = fid.readline().rstrip()

    major_version = line[0:9].strip().split('.')[0]

    if major_version == '2':
        result = readRinexObs211(
            filename, readSS=None, readLLI=None,
            includeAllGNSSsystems=None, includeAllObsCodes=None,
            desiredGNSSsystems=desiredGNSSsystems,
            desiredObsCodes=None, desiredObsBands=None,
        )
        v3_path_handled_decimation = False
    else:
        result = readRinexObs304(
            filename, readSS, readLLI,
            includeAllGNSSsystems, includeAllObsCodes,
            desiredGNSSsystems, desiredObsCodes, desiredObsBands,
            desired_data_rate=desired_data_rate,
        )
        v3_path_handled_decimation = True

    if isinstance(result, RinexObsData):
        data = result
    else:
        # Legacy reader still returns a raw tuple — wrap it
        data = RinexObsData.from_tuple(result)

    # v2 path: apply post-decimation fallback (parser-level skip not supported there)
    if desired_data_rate is not None and not v3_path_handled_decimation:
        data = _decimate_rinex_obs_data(data, float(desired_data_rate))
    elif desired_data_rate is not None and desired_data_rate <= 0:
        # Validate even when v3 path didn't get to (e.g. tInterval was NaN)
        raise ValueError(
            f"desired_data_rate must be positive (got {desired_data_rate})"
        )
    return data


def _decimate_rinex_obs_data(data: 'RinexObsData', desired_data_rate: float) -> 'RinexObsData':
    """Down-sample :class:`RinexObsData` to a coarser data rate.

    Keeps every ``stride``-th epoch where ``stride = round(desired_data_rate
    / tInterval)``.  When the desired rate is finer (smaller) than the
    file's native rate, the data is returned unchanged.

    Parameters
    ----------
    data : RinexObsData
        Parsed RINEX observation data.
    desired_data_rate : float
        Desired interval in seconds (must be > 0).

    Returns
    -------
    RinexObsData
        New object containing only the retained epochs.  ``tInterval``,
        ``nepochs`` and ``tLastObs`` are updated to reflect the
        decimated data.
    """
    if desired_data_rate <= 0:
        raise ValueError(
            f"desired_data_rate must be positive (got {desired_data_rate})"
        )

    t_interval = data.tInterval
    if t_interval is None or (isinstance(t_interval, float) and np.isnan(t_interval)) or t_interval <= 0:
        # Cannot decimate without a known native interval
        return data

    stride = int(round(desired_data_rate / t_interval))
    if stride <= 1:
        # Desired rate is finer than (or equal to) native rate — nothing to do
        return data

    nepochs = int(data.nepochs)
    if nepochs == 0:
        return data

    # 0-based indices of epochs to keep
    keep_idx = list(range(0, nepochs, stride))
    if not keep_idx:
        return data

    # ── Re-key per-system epoch dicts (1-based) ──
    def _reindex(epoch_dict_per_sys):
        if not isinstance(epoch_dict_per_sys, dict):
            return epoch_dict_per_sys
        new = {}
        for sys, epoch_dict in epoch_dict_per_sys.items():
            if not isinstance(epoch_dict, dict):
                new[sys] = epoch_dict
                continue
            new_sys = {}
            for new_e_idx, old_idx in enumerate(keep_idx, start=1):
                old_key = old_idx + 1  # 1-based original epoch key
                if old_key in epoch_dict:
                    new_sys[new_e_idx] = epoch_dict[old_key]
            new[sys] = new_sys
        return new

    new_GNSS_obs = _reindex(data.GNSS_obs)
    new_GNSS_LLI = _reindex(data.GNSS_LLI)
    new_GNSS_SS  = _reindex(data.GNSS_SS)

    # ── Decimate GNSS_SVs (per-system 2D ndarray, rows = epochs) ──
    new_GNSS_SVs = {}
    if isinstance(data.GNSS_SVs, dict):
        for sys, arr in data.GNSS_SVs.items():
            if isinstance(arr, np.ndarray) and arr.ndim == 2 and arr.shape[0] == nepochs:
                new_GNSS_SVs[sys] = arr[keep_idx, :].copy()
            else:
                new_GNSS_SVs[sys] = arr
    else:
        new_GNSS_SVs = data.GNSS_SVs

    # ── Decimate time_epochs (rows = epochs) ──
    if isinstance(data.time_epochs, np.ndarray) and data.time_epochs.ndim == 2 and data.time_epochs.shape[0] == nepochs:
        new_time_epochs = data.time_epochs[keep_idx, :].copy()
    else:
        new_time_epochs = data.time_epochs

    # ── Compute new metadata ──
    new_nepochs = len(keep_idx)
    new_tInterval = float(t_interval * stride)

    # tLastObs: derive from last kept epoch in time_epochs if possible
    new_tLastObs = data.tLastObs
    if isinstance(new_time_epochs, np.ndarray) and new_time_epochs.size and new_time_epochs.shape[0] >= 1:
        try:
            from gnssmultipath.Geodetic_functions import gpstime2date
            last_week = float(new_time_epochs[-1, 0])
            last_tow  = float(new_time_epochs[-1, 1])
            ymd_hms = gpstime2date(last_week, last_tow)  # [Y, M, D, H, M, S]
            new_tLastObs = np.array([[float(v)] for v in ymd_hms])
        except Exception:
            pass

    return RinexObsData(
        GNSS_obs=new_GNSS_obs,
        GNSS_LLI=new_GNSS_LLI,
        GNSS_SS=new_GNSS_SS,
        GNSS_SVs=new_GNSS_SVs,
        time_epochs=new_time_epochs,
        nepochs=new_nepochs,
        GNSSsystems=data.GNSSsystems,
        obsCodes=data.obsCodes,
        approxPosition=data.approxPosition,
        max_sat=data.max_sat,
        tInterval=new_tInterval,
        markerName=data.markerName,
        rinexVersion=data.rinexVersion,
        recType=data.recType,
        timeSystem=data.timeSystem,
        leapSec=data.leapSec,
        gnssType=data.gnssType,
        rinexProgr=data.rinexProgr,
        rinexDate=data.rinexDate,
        antDelta=data.antDelta,
        tFirstObs=data.tFirstObs,
        tLastObs=new_tLastObs,
        clockOffsetsON=data.clockOffsetsON,
        GLO_Slot2ChannelMap=data.GLO_Slot2ChannelMap,
        success=data.success,
    )


def readRinexObs304(filename, readSS=None, readLLI=None, includeAllGNSSsystems=None,includeAllObsCodes=None, \
                    desiredGNSSsystems=None, desiredObsCodes=None, desiredObsBands=None, desired_data_rate=None):
    """
    Program/function to read GNSS observations in RINEX 3.04/4.xx observation files.

    RINEX v4 observation files use the same record format as v3.xx,
    so this reader handles both versions.

    The main core of the program is 4 functions:
                                  rinexReadObsFileHeader304
                                  rinexReadObsBlockHead304
                                  rinexReadObsBlock304
                                  rinexFindNEpochs304


    To export every parameter use this code:

    GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
         obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
         rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success = \
         readRinexObs304(rinObsFilename)

    Tip: Unwanted outputs can be ignored using _

    --------------------------------------------------------------------------------------------------------------------------
    INPUTS

    filename:                 path and name of RINEX 3.04 observation file,
                              string

    readSS:                   Boolean, 0 or 1.
                              1 = read "Signal Strength" Indicators
                              0 = do not read "Signal Strength" Indicators

    readLLI:                  Boolean, 0 or 1.
                              1 = read "Loss-Of-Lock Indicators"
                              0 = do not read "Loss-Of-Lock Indicators"

    includeAllGNSSsystems:    Boolean, 0 or 1.
                              1 = include alle GNSS systems(GPS, GLONASS, Galieo, BeiDou)
                              0 = include only GNSSsystems specified in desiredGNSSsystems

    includeAllObsTypes:       Boolean, 0 or 1.
                              1 = include all valid ObsTypes
                              0 = include only ObsTypes specified in desiredObsTypes

    desiredGNSSsystems:       array og strings containing  codes of desired
                              GNSSsystems to be included,
                              ex. ["G", "E"]
                              OBS: Must be string array, NOT char vector

    desiredObsTypes:          array of strings containing desired ObsTypes to be
                              included, ex. ["C", "L", "S", "D"]
                              OBS: Must be string array, NOT char vector

    desiredObsBands:          array of desired obs Bands to be included,
                              ex [1, 5]
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS

    GNSS_obs:                 cell containing a matrix for each GNSS system.
                              Each matrix is a 3D matrix containing all
                              observation of current GNSS system for all epochs.
                              Order of obsType index is same order as in
                              obsCodes cell

                              GNSS_obs{GNSSsystemIndex}(PRN, obsType, epoch)
                                              GNSSsystemIndex: double,
                                              1,2,...,numGNSSsystems
                                              PRN: double
                                              ObsType: double: 1,2,...,numObsTypes
                                              epoch: double

    GNSS_LLI:                 cell containing a matrix for each GNSS system
                              Each matrix stores loss of lock indicators for
                              each epoch for that GNSS system

    GNSS_SS:                  cell containing a matrix for each GNSS system
                              Each matrix stores signal strength indicators
                              for each epoch for that GNSS system

    GNSS_SVs:                 cell containing a matrix for each GNSS system.
                              Each matrix contains number of satellites with
                              obsevations for each epoch, and PRN for those
                              satellites

                              GNSS_SVs{GNSSsystemIndex}(epoch, j)
                                              j=1: number of observed satellites
                                              j>1: PRN of observed satellites

    time_epochs:              matrix conatining gps-week and "time of week"
                              for each epoch
                              time_epochs(epoch,i),   i=1: week
                                                      i=2: time-of-week in seconds (tow)

    nepochs:                  number of epochs with observations in rinex observation file.

    GNSSsystems:              cell array containing codes of GNSS systems included
                              in RINEX observationfile. Elements are strings.
                              ex. "G" or "E"

    obsCodes:                 Cell that defines the observation
                              codes available for all GNSS system. Each cell
                              element is another cell containing the codes for
                              that GNSS system. Each element in this cell is a
                              string with three-characters. The first
                              character (a capital letter) is an observation code
                              ex. "L" or "C". The second character (a digit)
                              is a frequency code. The third character(a Capital letter)
                              is the attribute, ex. "P" or "X"

    approxPosition:           array containing approximate position from rinex
                              observation file header. [X, Y, Z]

    max_sat:                  array conataining max PRN number for each GNSS
                              system. Follows same order as GNSSsystems

    tInterval:                observations interval; seconds.

    markerName:               name of the antenna marker; '' if not specified

    rinexVersion:             string. rinex observation file version

    recType:                  receiver type, char vector

    timeSystem:               three-character code string of the time system
                              used for expressing tfirstObs;
                              can be GPS, GLO or GAL;

    leapSec:                  number of leap seconds since 6-Jan-1980.
                              UTC=GPST-leapSec. NaN by default.
                              THIS IS RINEX 3.04 OPTIONAL DATA

                              rinexHeader: cell column-vector containing the
                              following data:
                              rinexVersion:   RINEX version number; string.
                                              '' if not specified
                              rinexType:      RINEX file type; char

    gnssType:                 GNSS system of the satellites observed; can be
                              'G', 'R', 'E', 'C' or 'M' that stand for
                              GPS, GLONASS, GALILEO, BeiDou or Mixed ; char

    rinexProgr:               name of the software used to produce de RINEX
                              GPS obs file; '' if not specified


    rinexDate:                date/time of the RINEX file creation; '' if not
                              specified; char


    antDelta:                 column vector ot the three components of the
                              distance from the marker to the antenna,
                              in the following order - up, east and north;
                              null vector by default

    tFirstObs:                time stamp of the first observation record in the RINEX
                              observations file; column vector
                              [YYYY; MM; DD; hh; mm; ss.sssssss];
                              THIS IS CRITICAL DATA

    tLastObs:                 time stamp of the last observation record in the RINEX
                              observations file; column vector
                              [YYYY, MM, DD, hh, mm,ss.sssssss]. NaN by default.
                              THIS IS RINEX 3.04 OPTIONAL DATA

    clockOffsetsON:           receiver clock offsets flag. O if no realtime-derived
                              receiver clock offset was applied to epoch,
                              code and phase data (in other words, if the
                              file only has raw data), 1 otherwise. 0 by default.
                              THIS IS RINEX 3.04 OPTIONAL DATA

    GLO_Slot2ChannelMap:      map container that maps GLONASS slot numbers to
                              their respective channel number.
                              GLO_Slot2ChannelMap(slotnumber)

    success:                  Boolean. 1 if the reading of the RINEX
                              observations file seems to be successful,
                              0 otherwise
    --------------------------------------------------------------------------------------------------------------------------

    ADVICE: The function rinexFindNEpochs() calculates the amount of observation epochs in
    advance. This calculation will be incredibly more effective if TIME OF
    LAST OBS is included in the header of the observation file. It is
    strongly advized to manually add this information to the header if it is
    not included by default.
    --------------------------------------------------------------------------------------------------------------------------

    According to RINEX 3.04 the observation type codes are:
    Observation code
            C: Pseudorange
    GPS: C/A, L2C
    Glonass: C/A
    Galileo: All
    L: Carrier phase
    D: Doppler frequency
    S: Raw signal strengths or SNR values as given by the receiver for the
    respective phase observations
    I: Ionosphere phase delay
    X: Receiver channel numbers

    Frequency code
    GPS Glonass Galileo SBAS
    1: L1 G1 E1 B1    (GPS,QZSS,SBAS,BDS)
    2: L2 G2 B1-2     (GLONASS)
    4: G1a            (Galileo)
    5: L5 E5a B2/B2a  (GPS, QZSS, SBAS, IRNSS)
    6: L6 E6 B3 G2a   (Galileo, QZSS, BDS, GLONASS)
    7: E5b B2/B2b     (Galileo)
    8: E5a+b E5a+b    (Galileo, BDS)
    9: S              (IRNSS)
    0: for type X     (all)

    Attribute:
    A = A channel     (Galileo,IRNSS,GLONASS)
    B = B channel     (Galileo,IRNSS,GLONASS)
    C = C channel     (Galiloe, IRNSS)
    C code-based  (SBAS, GPS, GLONASS, QZSS)
    D = Semi-codelss  (GPS)

    I = I channel     (GPS, Galileo, QZSS, BDS)
    L = L channel     (L2C GPS, QZSS)
    P channel     (GPS. QZSS)
    M = M code-based  (GPS)
    N = Codeless      (GPS)
    P = P code-based  (GPS, GLONASS)
    Pilot channel (BDS)

    Q = Q channel     (GPS, Galileo, QZSS, BDS)
    S = D channel     (GPS, Galileo, QZSS, BDS)
    M channel     (L2C GPS, QZSS)

    W = Based on Z-tracking (GPS)
    X = B+C channels  (Galileo, IRNSS)
    I+Q channels  (GPS, IRNSS)
    M+L channels  (GPS, QZSS)
    D+P channels  (GPS, QZSS, BDS)

    Y = Y code based  (GPS)
    Z = A+B+C channels(Galileo)
    D+P channels  (BDS)
    --------------------------------------------------------------------------------------------------------------------------
    """
    ## -- Setting None arguments
    if readSS is None:
        readSS = 1
    if readLLI is None:
        readLLI = 1
    if includeAllGNSSsystems is None and desiredGNSSsystems is None:
        includeAllGNSSsystems = 1
    if includeAllObsCodes is None and desiredObsCodes is None:
        includeAllObsCodes = 1
    if desiredGNSSsystems is None:
        desiredGNSSsystems = ['G','R','E','C']
    if desiredObsCodes is None:
        desiredObsCodes = ['C','L','S','D']
    if desiredObsBands is None:
        desiredObsBands = list(np.arange(1,10))

    ## Get the start time
    t = time.process_time()

    ## - Initialize variables in case of input error
    GNSS_obs       = np.nan
    GNSS_LLI       = np.nan
    GNSS_SS        = np.nan
    GNSS_SVs       = np.nan
    time_epochs    = np.nan
    nepochs        = np.nan
    GNSSsystems    = np.nan
    obsCodes       = np.nan
    approxPosition = np.nan
    max_sat        = np.nan
    tInterval      = np.nan
    markerName     = np.nan
    rinexVersion   = np.nan
    recType        = np.nan
    timeSystem     = np.nan
    leapSec        = np.nan
    gnssType       = np.nan
    rinexProgr     = np.nan
    rinexDate      = np.nan
    antDelta       = np.nan
    tFirstObs      = np.nan
    tLastObs       = np.nan
    clockOffsetsON = np.nan
    GLO_Slot2ChannelMap = np.nan

    ### --- Dict for storing data
    GNSS_obs = {}
    GPS = {}
    GLONASS = {}
    Galileo = {}
    BeiDou = {}
    GPS_LLI = {}
    GLONASS_LLI = {}
    Galileo_LLI = {}
    BeiDou_LLI = {}
    GPS_SS = {}
    GLONASS_SS = {}
    Galileo_SS = {}
    BeiDou_SS = {}

    ## -- Test if readSS is boolean
    if readSS!=1 and readSS!=0:
        print('INPUT ERROR(readRinexObs304): The input argument readSS must be either 1 or 0')
        success = 0
        return


    ## -- Test if readLLI is boolean
    if readLLI!=1 and readLLI!=0:
        print('INPUT ERROR(readRinexObs304): The input argument readLLI must be either 1 or 0')
        success = 0
        return

    max_GPS_PRN     = 36 # Max number of GPS PRN in constellation
    max_GLONASS_PRN = 36 # Max number of GLONASS PRN in constellation
    max_Galileo_PRN = 36 # Max number of Galileo PRN in constellation
    max_Beidou_PRN  = 100 # Max number of BeiDou PRN in constellation

    ## -- Read header of observation file
    [success, rinexVersion, gnssType, markerName, recType, antDelta,\
    GNSSsystems,numOfObsCodes, obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval, \
    timeSystem, _, clockOffsetsON, rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, _, fid] = \
    rinexReadObsFileHeader304(filename, includeAllGNSSsystems, includeAllObsCodes,desiredGNSSsystems, desiredObsCodes, desiredObsBands)

    if success==0:
        return

    ## -- Read all remaining observation data at once for fast processing
    remaining_lines = fid.readlines()
    fid.close()

    ## -- Count epochs and determine timing from buffered lines
    epoch_line_indices = [i for i, ln in enumerate(remaining_lines) if ln.startswith('>')]
    nepochs = len(epoch_line_indices)

    if nepochs == 0:
        print('ERROR(readRinexObs304): No epoch records found in observation file')
        success = 0
        return

    # Compute tInterval from first two epoch headers if not in header
    if np.isnan(tInterval) and nepochs >= 2:
        _t1 = [float(x) for x in remaining_lines[epoch_line_indices[0]][1:].split() if x][:6]
        _t2 = [float(x) for x in remaining_lines[epoch_line_indices[1]][1:].split() if x][:6]
        tInterval = _t2[5] - _t1[5]

    # Compute tLastObs if not in header
    try:
        _tl_is_nan = np.all(np.isnan(tLastObs))
    except (TypeError, ValueError):
        _tl_is_nan = False
    if _tl_is_nan:
        _lp = [float(x) for x in remaining_lines[epoch_line_indices[-1]][1:].split() if x][:6]
        tLastObs = np.array([[int(_lp[0])], [int(_lp[1])], [int(_lp[2])],
                            [int(_lp[3])], [int(_lp[4])], [int(_lp[5])]])
        print('INFO(readRinexObs304): TIME OF LAST OBS computed from observation data')

    ## -- Determine decimation stride for parser-level down-sampling
    decim_stride = 1
    if desired_data_rate is not None:
        if desired_data_rate <= 0:
            raise ValueError(
                f"desired_data_rate must be positive (got {desired_data_rate})"
            )
        if not (isinstance(tInterval, float) and np.isnan(tInterval)) and tInterval > 0:
            s = int(round(float(desired_data_rate) / float(tInterval)))
            if s > 1:
                decim_stride = s
                # nepochs becomes the count of epochs we will actually keep
                nepochs = (nepochs + decim_stride - 1) // decim_stride
                tInterval = float(tInterval * decim_stride)
                print(
                    f'INFO(readRinexObs304): Decimating observations by stride {decim_stride} '
                    f'(desired_data_rate={desired_data_rate}s -> tInterval={tInterval}s, '
                    f'kept {nepochs} epochs)'
                )

    success = 1
    nGNSSsystems = len(GNSSsystems)
    GNSS_SVs = {}
    max_sat  =  np.zeros([nGNSSsystems,1])
    t_week = []
    t_tow = []
    GNSSsystems_full_names =  [""]*nGNSSsystems
    GNSS_LLI = {}
    GNSS_SS = {}

    # Create array for max_sat. Initialize cell elements in dicts
    for k in np.arange(0,nGNSSsystems):
        if GNSSsystems[k+1] == 'G':
            max_sat[k] = max_GPS_PRN
            GNSS_SVs['G'] = np.zeros([nepochs, int(max_sat[k].item() + 1)], dtype=np.int16)
            GNSSsystems_full_names[k] = "GPS"
        elif GNSSsystems[k+1] == 'R':
            max_sat[k] = max_GLONASS_PRN
            GNSS_SVs['R'] = np.zeros([nepochs, int(max_sat[k].item() + 1)], dtype=np.int16)
            GNSSsystems_full_names[k] = "GLONASS"

        elif GNSSsystems[k+1] == 'E':
            max_sat[k] = max_Galileo_PRN
            GNSS_SVs['E'] = np.zeros([nepochs, int(max_sat[k].item() + 1)], dtype=np.int16)
            GNSSsystems_full_names[k] = "Galileo"

        elif GNSSsystems[k+1] == 'C':
            max_sat[k] = max_Beidou_PRN
            GNSS_SVs['C'] = np.zeros([nepochs, int(max_sat[k].item() + 1)], dtype=np.int16)
            GNSSsystems_full_names[k] = "BeiDou"
        else:
            print(f'ERROR(readRinexObs304): Only following GNSS systems are compatible with this program: GPS, GLONASS, Galileo, Beidou. {GNSSsystems[k]} is not valid')
            return


        curr_sys = GNSSsystems[k+1]
        # Note: GNSS_obs[curr_sys] is populated as a per-epoch dict below
        # (e.g. GNSS_obs['G'] = GPS), so no large 3D pre-allocation is needed here.


        # Preallocation LLI and SS
        if readLLI:
            GNSS_LLI[curr_sys] = np.zeros([nGNSSsystems,1])
        else:
            GNSS_LLI[curr_sys] = np.nan
        if readSS:
            GNSS_SS[curr_sys] = np.zeros([nGNSSsystems,1])
        else:
            GNSS_SS[curr_sys] = np.nan


    GNSS_names = dict(zip(['G', 'R', 'E', 'C'],['GPS','GLONASS','Galileo','Beidou']))
    current_epoch      = 0
    file_epoch_index   = -1  # 0-based count of epochs encountered in file (incl. skipped)

    ## -- Pre-compute system-to-index lookup (avoids list comprehension per satellite)
    sys_to_idx = {v: k for k, v in GNSSsystems.items()}

    ## -- Pre-compute max_sat as plain ints and obs code field positions per system
    max_sat_int = {k+1: int(max_sat[k].item()) if hasattr(max_sat[k], "item") else int(max_sat[k]) for k in range(nGNSSsystems)}
    sys_obs_field_positions = {}
    for k in range(1, nGNSSsystems + 1):
        sys_obs_field_positions[k] = [4 + ci * 16 for ci in obsCodeIndex[k]]

    ## -- Per-system storage dicts (mapped by system char for fast access)
    _obs_dicts = {'G': GPS, 'R': GLONASS, 'E': Galileo, 'C': BeiDou}
    _lli_dicts = {'G': GPS_LLI, 'R': GLONASS_LLI, 'E': Galileo_LLI, 'C': BeiDou_LLI}
    _ss_dicts  = {'G': GPS_SS, 'R': GLONASS_SS, 'E': Galileo_SS, 'C': BeiDou_SS}

    ## -- Initialize progress bar
    n_update_break = max(int(nepochs // 10), 1)
    bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'
    n_remaining = len(remaining_lines)

    with tqdm(total=100, desc="Rinex observations are being read", position=0, leave=True, bar_format=bar_format) as pbar:
        line_cursor = 0
        while line_cursor < n_remaining:
            line = remaining_lines[line_cursor]
            line_cursor += 1

            ## -- Find next epoch header line (starts with '>')
            if not line.startswith('>'):
                continue
            line = line.rstrip()

            ## -- Handle event flags (epochflag > 1): skip special records
            epochflag = int(line[31])
            while epochflag > 1:
                linejump = int(line[32:35])
                for _ in range(linejump + 1):
                    if line_cursor >= n_remaining:
                        break
                    line = remaining_lines[line_cursor].rstrip()
                    line_cursor += 1
                if line_cursor >= n_remaining:
                    break
                epochflag = int(line[31])

            if line_cursor >= n_remaining and epochflag > 1:
                break

            numSV = int(line[32:35])

            ## -- Parser-level decimation: skip non-kept epochs entirely
            file_epoch_index += 1
            if decim_stride > 1 and (file_epoch_index % decim_stride) != 0:
                # Advance past this epoch's satellite lines without parsing
                line_cursor += numSV
                continue

            ## -- Parse epoch date
            date = [float(el) for el in line[1:].split() if el][:6]

            current_epoch += 1
            if current_epoch % n_update_break == 0:
                pbar.update(10)

            ## Convert date to GPS-week and "time-of-week"
            week, tow = date2gpstime(int(date[0]), int(date[1]), int(date[2]), int(date[3]), int(date[4]), int(date[5]))
            t_week.append(week)
            t_tow.append(tow)

            ## -- Initialize per-epoch arrays
            nGNSS_sat_current_epoch = [0] * nGNSSsystems
            GNSS_obs_dum = {}
            GNSS_LLI_dum = {}
            GNSS_SS_dum  = {}
            for k in range(nGNSSsystems):
                mk = max_sat_int[k+1]
                nc = numOfObsCodes[k]
                GNSS_obs_dum[k+1] = np.zeros((mk + 1, nc))
                GNSS_LLI_dum[k+1] = np.zeros((mk + 1, nc))
                GNSS_SS_dum[k+1]  = np.zeros((mk + 1, nc))

            ## -- Parse satellite observation lines for this epoch
            for sat in range(numSV):
                if line_cursor >= n_remaining:
                    break
                sat_line = remaining_lines[line_cursor].rstrip()
                line_cursor += 1
                if not sat_line:
                    continue

                sys_char = sat_line[0]
                if sys_char not in sys_to_idx:
                    continue

                gi = sys_to_idx[sys_char]  # GNSSsystem index (1-based)
                prn = int(sat_line[1:3])
                nGNSS_sat_current_epoch[gi - 1] += 1

                n_obs = numOfObsCodes[gi - 1]
                field_positions = sys_obs_field_positions[gi]

                # Pad line to ensure safe field extraction
                max_pos_needed = field_positions[-1] + 16 if field_positions else 4
                padded = sat_line.ljust(max_pos_needed)

                # Extract all observation values, LLI, and SS from the line
                obs_row = GNSS_obs_dum[gi][prn]
                lli_row = GNSS_LLI_dum[gi][prn] if readLLI else None
                ss_row  = GNSS_SS_dum[gi][prn] if readSS else None

                for obs_num in range(n_obs):
                    pos = field_positions[obs_num]
                    val_str = padded[pos:pos+14].strip()
                    obs_row[obs_num] = float(val_str) if val_str else 0.0

                    if readLLI:
                        lli_ch = padded[pos+13]
                        lli_row[obs_num] = int(lli_ch) if not lli_ch.isspace() else -999

                    if readSS:
                        ss_ch = padded[pos+14]
                        ss_row[obs_num] = int(ss_ch) if not ss_ch.isspace() else -999

                # Store PRN in GNSS_SVs
                GNSS_SVs[sys_char][current_epoch-1, nGNSS_sat_current_epoch[gi-1]] = prn

            ## -- Store per-system data for this epoch
            for k in range(nGNSSsystems):
                curr_sys = GNSSsystems[k+1]
                GNSS_SVs[curr_sys][current_epoch-1, 0] = nGNSS_sat_current_epoch[k]

                _obs_dicts[curr_sys][current_epoch] = GNSS_obs_dum[k+1]
                if readLLI:
                    _lli_dicts[curr_sys][current_epoch] = GNSS_LLI_dum[k+1]
                if readSS:
                    _ss_dicts[curr_sys][current_epoch] = GNSS_SS_dum[k+1]

    ## -- Build time_epochs array (once, at end)
    time_epochs = np.column_stack((t_week, t_tow)) if t_week else np.empty((0, 2))

    ## -- When parser-level decimation is active, recompute tLastObs
    ##    from the last kept epoch (header value referred to original last epoch)
    if decim_stride > 1 and time_epochs.shape[0] > 0:
        try:
            from gnssmultipath.Geodetic_functions import gpstime2date
            last_week = float(time_epochs[-1, 0])
            last_tow  = float(time_epochs[-1, 1])
            ymd_hms = gpstime2date(last_week, last_tow)
            tLastObs = np.array([[float(v)] for v in ymd_hms])
        except Exception:
            pass

    ## -- Storing observation in dictionary
    GNSS_obs['G'] = GPS
    GNSS_obs['R'] = GLONASS
    GNSS_obs['E'] = Galileo
    GNSS_obs['C'] = BeiDou
    ## -- Storing "loss of lock indicaors"  in dict
    GNSS_LLI['G'] = GPS_LLI
    GNSS_LLI['R'] = GLONASS_LLI
    GNSS_LLI['E'] = Galileo_LLI
    GNSS_LLI['C'] = BeiDou_LLI
    ## -- Storing  SS in dict
    GNSS_SS['G'] = GPS_SS
    GNSS_SS['R'] = GLONASS_SS
    GNSS_SS['E'] = Galileo_SS
    GNSS_SS['C'] = BeiDou_SS
    # Deleting systems with no observations
    del_sys = list(GNSS_obs.keys())
    for sys in del_sys:
        if not GNSS_obs[sys]:
            del GNSS_obs[sys]

    if current_epoch!= nepochs and success == 1:
        print('ERROR(readRinexObs304): The amount of epochs calculated in advance(nepochs = %d) does not equal number og epochs prossesed(current_epoch = %d).\nCheck that header information concerning TIME OF FIRST OBS and TIME OF LAST OBS is correct.\n' %(nepochs, current_epoch))

    messages = {}
    if success == 1:
        messages[0]= 'INFO(readRinexObs304): The following GNSS systems have been read into the data:'
        for k in np.arange(0,nGNSSsystems):
            messages[k+1]= 'INFO(readRinexObs304): The following %s observation types have been registered:' % (GNSS_names[GNSSsystems[k+1]])
            curr_sys = GNSSsystems[k+1]
            for obs in np.arange(0, len(obsCodes[k+1][curr_sys])):
                if obs == 0:
                    messages[k+1]= messages[k+1] + ' %s' % (obsCodes[k+1][curr_sys][obs])
                else:
                    messages[k+1]= messages[k+1] + ', %s' % (obsCodes[k+1][curr_sys][obs])
            if k == 0:
                messages[0]= messages[0] + ' %s' % GNSS_names[GNSSsystems[1]]
            else:
                messages[0]= messages[0] + ', %s' % GNSS_names[GNSSsystems[k+1]]

        for msg in np.arange(0,len(messages)):
            print(messages[msg])

    if readLLI:
        print('INFO(readRinexObs304): LLI have been read (if present in observation file)')
    else:
        print('INFO(readRinexObs304): LLI have not been read')


    if readSS:
        print('INFO(readRinexObs304): SS have been read (if present in observation file)')
    else:
        print('INFO(readRinexObs304): SS have not been read')


    ## --  Finding processing time
    et = time.process_time()  # get the end time
    e = et - t                # get execution time

    if e >= 3600:
        hours = np.floor(e/3600)
        minutes = np.floor((e-hours*3600)/60)
        seconds = e-hours*3600-minutes*60
        print('INFO(readRinexObs304): Total processing time: %d hours, %d minutes, %f seconds\n' % (hours, minutes, seconds))
    elif e>60:
        minutes = np.floor(e/60)
        seconds = e-minutes*60
        print('INFO(readRinexObs304): Total processing time: %d minutes, %f seconds\n' % (minutes, seconds))
    else:
        print('INFO(readRinexObs304): Total processing time: %f seconds\n\n' % (e))


    return GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
        obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
        rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success


def rinexFindNEpochs304(filename, tFirstObs, tLastObs, tInterval):
    """
    Function that computes number of epochs in Rinex 3.xx observation file.
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS

    filename:         RINEX observation filename

    tFirstObs:        time stamp of the first observation record in the RINEX
                      observations file; column vector
                      [YYYY; MM; DD; hh; mm; ss.sssssss];

    tLastObs:         time stamp of the last observation record in the RINEX
                      observations file; column vector
                      [YYYY; MM; DD; hh; mm; ss.sssssss]. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function

    tInterval:        observations interval; seconds. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function.
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    -------

    nepochs:          number of epochs in Rinex observation file with
                      observations

    tLastObs:         time stamp of the last observation record in the RINEX
                      observations file; column vector
                      [YYYY, MM, DD, hh, mm, ss.sssssss]. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function

    tInterval:        observations interval; seconds. If this information
                      was not available in rinex observation header the
                      default value is Nan.

    success:                  Boolean. 1 if the function seems to be successful,
                              0 otherwise
    --------------------------------------------------------------------------------------------------------------------------

    ADVICE: The function rinexFindNEpochs() calculates the amount of observation epochs in
    advance. This calculation will be incredibly more effective if TIME OF
    LAST OBS is included in the header of the observation file. It is
    strongly advized to manually add this information to the header if it is
    not included by default.
    --------------------------------------------------------------------------------------------------------------------------
    """
    # #  Testing input arguments
    success = 1
    nepochs = 0

    ## --Test if filename is valid format
    if type(filename) is not str:
        raise TypeError('INPUT ERROR(rinexFindNEpoch): The input argument filename'\
            'is of type %s. Must be of type string' % type(filename))


    ## --Open observation file
    fid = open(filename, 'rt')
    #  tLastObs is in header
    if ~np.all(np.isnan(tLastObs)):
        #  tInterval is not in header
        if np.isnan(tInterval):
            tInterval_found = 0
            first_epoch_found = 0
            #  calculate tInterval
            while not tInterval_found: #  calculate tInterval
                line = fid.readline().rstrip()
                #  start of new epoch
                if '>' in line:
                    if not first_epoch_found: #  first epoch
                        first_epoch_time = line[1::]
                        first_epoch_time = [float(el) for el in line[1::].split(" ") if el != ""]
                        first_epoch_time = first_epoch_time[:6]
                        first_epoch_found = 1
                    else: #  seconds epoch
                        second_epoch_time = line[1::]
                        second_epoch_time = [float(el) for el in line[1::].split(" ") if el != ""]
                        second_epoch_time = second_epoch_time[:6]
                        tInterval = second_epoch_time[5]-first_epoch_time[5]
                        tInterval_found = 1

        fid.close(); fid = open(filename, 'rt')
        tFirstObs = tFirstObs.astype(int)
        tLastObs = tLastObs.astype(int)
        rinex_lines = fid.readlines()
        epoch_line = [line for line in rinex_lines if line.startswith('>')] # list with all the line thats defines a epoch
        nepochs = len(epoch_line)
    #  if tLastObs is not in header. Function counts number of epochs manually
    else:
    ## New code for finding last
        print('INFO(rinexFindEpochs304): The header of the rinex observation file does not contain TIME OF LAST OBS.\n' \
            'This will be calculated, but consider editing rinex header to include TIME OF LAST HEADER')

        fid.close(); fid = open(filename, 'rt')
        rinex_lines = fid.readlines()
        epoch_lines = [line for line in rinex_lines if '>' in line] # list with all the line thats defines a epoch
        nepochs = len(epoch_lines)
        ## Computing the tInterval if not present in the header
        if np.isnan(tInterval):
            first_epoch_line  = epoch_lines[0][1::]
            second_epoch_line = epoch_lines[1][1::]
            first_epoch_time  = [float(el) for el in first_epoch_line[1::].split(" ") if el != ""]
            second_epoch_time = [float(el) for el in second_epoch_line[1::].split(" ") if el != ""]
            tInterval = second_epoch_time[5]-first_epoch_time[5]

        line = epoch_lines[-1]
        line = line[1:60]     #  deletes 'TIME OF LAST OBS'
        line_ = [el for el in line.split(" ") if el != ""]
        for k in np.arange(0,6):
            tok = line_.pop(0)
            if k ==0:
                yyyy = int(tok)
            elif k ==1:
                mm = int(tok)
            elif k ==2:
                dd = int(tok)
            elif k ==3:
                hh = int(tok)
            elif k ==4:
                mnt = int(tok)
            elif k ==5:
                ss = float(tok)

        tLastObs = np.array([[yyyy],[mm],[dd],[hh],[mnt],[ss]]).astype(int)
        print('INFO(rinexFindNEpochs304): TIME OF LAST OBS has been found and amount of epochs have been computed')
        fid.close()
    return int(nepochs), tLastObs, tInterval, success


def rinexReadObsFileHeader304(filename, includeAllGNSSsystems, includeAllObsCodes,desiredGNSSsystems,desiredObsCodes, desiredObsBands):
    """
    Extracts relevant data from the header of a RINEX 3.xx GNSS observations
    file. Excludes undesired GNSS systems, obsevation codes and/or frequency
    bands.

    --------------------------------------------------------------------------------------------------------------------------
    INPUTS:
    ------

    filename:                     RINEX observation filename and path

    includeAllGNSSsystems:        Boolean, 0 or 1.
                                      1 = include alle GNSS systems
                                          (GPS, GLONASS, Galieo, BeiDou)
                                      0 = include only GNSSsystems
                                          specified in desiredGNSSsystems

    includeAllobsCodes:           Boolean, 0 or 1.
                                      1 = include all valid obsCodes
                                      0 = include only obsCodes
                                          specified in desiredobsCodes

    desiredGNSSsystems:           string array containing desired GNSSsystems
                                  to be included, ex. ["G", "E", "C"]

    desiredobsCodes:              string array containing desired obsCodes to
                                  be included, ex. ["C", "L", "S", "D"]

    desiredObsBands:              array of desired obs Bands to be included,
                                  ex [1, 5]

    NOTE: If both includeAllGNSSsystems and includeAllobsCodes Boolean are 1
          then the last three input arguments are optional to include and may
          be left blank without en error.
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    -------

    success:                      1 if the reading of the RINEX observations
                                  file seems to be successful, 0 otherwise

    rinexVersion:                 string. rinex observation file version

    gnssType:                     GNSS system of the satellites observed; can
                                  be 'G', 'R', 'E', 'C' or 'M' that stand for
                                  GPS, GLONASS, GALILEO, BeiDou or Mixed; char

    markerName:                   name of the antenna marker; '' if not
                                  specified

    recType:                      Receiver type, char vector

    antDelta:                     column vector ot the three components of
                                  the distance from the marker to the antenna,
                                  in the following order - up, east and north;
                                  null vector by default

    GNSSsystems:                  cell array containing codes of GNSS systems
                                  included in RINEX observationfile. Elements
                                  are strings. ex. "G" or "E"

    numOfObsCodes:                column vector containing number of observation
                                  types for each GNSS system. Order is the same
                                  as GNSSsystems

    obsCodes:                     Cell that defines the observation
                                  codes available for all GNSS system. Each
                                  cell element is another cell containing the
                                  codes for that GNSS system. Each element in
                                  this cell is a string with three-characters.
                                  The first character (a capital letter) is
                                  an observation code ex. "L" or "C". The
                                  second character (a digit) is a frequency
                                  code. The third character(a Capital letter)
                                  is the attribute, ex. "P" or "X"

    obsCodeIndex:                 cell with one cell element for each GNSS
                                  system. Order is the same as GNSSsystems.
                                  Each cell element contains an array of
                                  indices. These indices indicate the
                                  observation types that should be read
                                  for each GNSS system. ex. If one index for
                                  GPS is 1 then the first observation type
                                  for GPS should  be read.

    tFirstObs:                    time stamp of the first observation record
                                  in the RINEX observations file; column vector
                                  [YYYY; MM; DD; hh; mm; ss.sssssss];
                                  THIS IS CRITICAL DATA

    tLastObs:                     time stamp of the last observation record
                                  in the RINEX observations file; column vector
                                  [YYYY; MM; DD; hh; mm;ss.sssssss].
                                  NaN by default.
                                  THIS IS RINEX 3.04 OPTIONAL DATA

    tInterval:                    observations interval; seconds.

    timeSystem:                   three-character code string of the time
                                  system used for expressing tfirstObs;
                                  can be GPS, GLO or GAL;

    numHeaderLines:               number of lines in header

    rinexProgr:                   name of the software used to produce de
                                  RINEX GPS obs file; '' if not specified

    rinexDate:                    date/time of the RINEX file creation; ''
                                  if not specified; char

    leapSec:                      number of leap seconds since 6-Jan-1980.
                                  UTC=GPST-leapSec. NaN by default.
                                  THIS IS RINEX 3.04 OPTIONAL DATA

    approxPosition:               array containing approximate position from
                                  rinex observation file header. [X, Y, Z]

    GLO_Slot2ChannelMap:          map container that maps GLONASS slot
                                  numbers to their respective channel number.
                                  GLO_Slot2ChannelMap(slotnumber)

    eof:                          end-of-file flag; 1 if end-of-file was reached,
                                  0 otherwise

    fid:                          Matlab file identifier of a Rinex
                                  observations text file
    --------------------------------------------------------------------------------------------------------------------------

    According to RINEX 3.04 these codes are:

       Observation codes:
       ------------------
       C: Pseudorange
          GPS: C/A, L2C
          Glonass: C/A
          Galileo: All
       L: Carrier phase
       D: Doppler frequency
       S: Raw signal strengths or SNR values as given by the receiver for the
          respective phase observations
       I: Ionosphere phase delay
       X: Receiver channel numbers

       Frequency code:
       ---------------
       GPS Glonass Galileo SBAS
       1: L1 G1 E1 B1    (GPS,QZSS,SBAS,BDS)
       2: L2 G2 B1-2     (GLONASS)
       4: G1a            (Galileo)
       5: L5 E5a B2/B2a  (GPS, QZSS, SBAS, IRNSS)
       6: L6 E6 B3 G2a   (Galileo, QZSS, BDS, GLONASS)
       7: E5b B2/B2b     (Galileo)
       8: E5a+b E5a+b    (Galileo, BDS)
       9: S              (IRNSS)
       0: for type X     (all)

       Attribute:
       ----------
       A = A channel     (Galileo,IRNSS,GLONASS)
       B = B channel     (Galileo,IRNSS,GLONASS)
       C = C channel     (Galiloe, IRNSS)
           C code-based  (SBAS, GPS, GLONASS, QZSS)
       D = Semi-codelss  (GPS)

       I = I channel     (GPS, Galileo, QZSS, BDS)
       L = L channel     (L2C GPS, QZSS)
           P channel     (GPS. QZSS)
       M = M code-based  (GPS)
       N = Codeless      (GPS)
       P = P code-based  (GPS, GLONASS)
           Pilot channel (BDS)

       Q = Q channel     (GPS, Galileo, QZSS, BDS)
       S = D channel     (GPS, Galileo, QZSS, BDS)
           M channel     (L2C GPS, QZSS)

       W = Based on Z-tracking (GPS)
       X = B+C channels  (Galileo, IRNSS)
           I+Q channels  (GPS, IRNSS)
           M+L channels  (GPS, QZSS)
           D+P channels  (GPS, QZSS, BDS)

       Y = Y code based  (GPS)
       Z = A+B+C channels(Galileo)
           D+P channels  (BDS)
    -------------------------------------------------------------------------------------------------------------------------
    """
    eof         = 0
    success     = 1
    warnings    = 0
    antDelta    = []
    timeSystem  = ''
    tFirstObs   = []
    tLastObs    = np.nan
    tInterval   = np.nan
    rinexProgr  = np.nan
    rinexDate   = np.nan
    obsCodes    = {}
    GNSSsystems = {}
    gnssType    = ""
    markerName  = ""
    numHeaderLines  = 0
    clockOffsetsON  = 0
    numGNSSsystems  = 0
    leapSec         = np.nan
    numOfObsCodes   = []
    rinexHeader     = {}
    approxPosition  = [0, 0, 0]
    obsCodeIndex = {}
    rinexVersion = np.nan
    recType = np.nan
    GLO_Slot2ChannelMap = np.nan

    ## -------Testing input arguments
    # Test if filename is valid format
    if type(filename) != str:
        print('INPUT ERROR(rinexReadsObsHeader304): The input argument filename is of type %s.\n Must be of type string or char' %(type(filename)))
        success = 0
        fid     = 0
        return success
    ## -- Open rinex observation file
    fid = open(filename,'r')
    if os.stat(filename).st_size == 0:
        raise ValueError('ERROR: This file seems to be empty')

    while 1: # Gobbling the header
        numHeaderLines = numHeaderLines + 1
        line = fid.readline().rstrip()
        if 'END OF HEADER' in line:
            break
        if numHeaderLines == 1: # if first line of header
            rinexVersion = line[0:9]
            # store rinex type, ex. "N" or "O"
            rinexType = line[20]
            # if rinex file is not an observation file
            if rinexType != 'O':  # Rinex file is oservation file
                print('ERROR(rinexReadObsFileHeader304): the file is not a RINEX observations data file!')
                success = 0
                fid.close()
                return

            ## -- Check gnss type  ## Changend indent here 09.12.2022 (was apart of the if test above earlier, and thats wrong)
            gnssType = line[40] # reads the GNSS system type
            if gnssType not in [' ', 'G', 'R', 'C', 'E', 'M' ]:
                if gnssType in ['J', 'I', 'S']:
                    print('ERROR(rinexReadObsFileHeader304): This software is meant for reading GNSS data only.\
                           %s is an invalid satellite system type.' %(gnssType))
                else:
                    print('ERROR(rinexReadObsFileHeader304): %s is an unrecognized satellite system type.' %(gnssType))

                success = 0
                fid.close()
            ## -- If no system type, set G
            if gnssType == ' ':
                gnssType = 'G'

        if 'PGM / RUN BY / DATE' in line:
            rinexProgr = line[0:20] # rinex program
            rinexDate = line[40:60] # rinex date

        if 'MARKER NAME' in line:
            markerName = line.strip() # markername

        ## if no marker name, "MARKER" is read, so set to blank
        if 'Marker' in markerName:
            markerName = ''

        if 'ANTENNA: DELTA H/E/N' in line:
            for k in np.arange(0,3):
                line_ = [el for el in line.split(" ") if el != ""]
                antDelta = [line_[0],line_[1],line_[2]]

        ## Section describing what GNSS systems are present, and their obs types
        if 'SYS / # / OBS TYPES' in line:
            line = line[0:60]     # deletes 'SYS / # / OBS TYPES'
            line_ = [el for el in line.split(" ") if el != ""]
            Sys = line_.pop(0) # assingning system to variable and removing it from the list
            if Sys not in ["G","R","E","C"]: # added this line 29.01.2023 to fix bug where Only one system and several lines with Obscodes in rinex file
                continue
            nObs = int(line_.pop(0))
            ## array for storing indeces of undesired ObsCodes for this GNSS system
            undesiredobsCodeIndex = []
            desiredObsCodeIndex = []
            ## is Sys amoung desired GNSS systems
            if (includeAllGNSSsystems and Sys in ["G", "R", "E", "C"] or Sys in desiredGNSSsystems):
                numGNSSsystems  = numGNSSsystems + 1 # increment number of GNSS systems
                GNSSsystems[numGNSSsystems] = str(Sys) # Store current GNSS system
                GNSSSystemObsCodes = {}  # Reset cell of obsCodes for this GNSS system
                obsCode_list = []
                for k in np.arange(0,nObs):
                    obsCode = line_.pop(0)
                    # Checking if obsCode is valid
                    if len(obsCode) != 3 or obsCode[0] not in ['C', 'L', 'D','S', 'I', 'X'] or  \
                              obsCode[1] not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'] or \
                               obsCode[2] not in ['A', 'B', 'C', 'D', 'I', 'L', 'M', 'N', 'P', 'Q', 'S', 'W', 'X', 'Y', 'Z']:
                        print('ERROR (rinexReadsObsHeader304):  obsCode %s is a not a standard RINEX 3.04 observation type!' %(obsCode))

                    ## is obsCode amoung desired obscodes and frequency bands
                    if includeAllObsCodes or obsCode[0] in desiredObsCodes and int(obsCode[1]) in desiredObsBands:
                         ## store obsCode if amoung desire obsCodes
                        obsCode_list.append(obsCode)
                        GNSSSystemObsCodes[Sys] =  obsCode_list
                        desiredObsCodeIndex.append(k)
                    else:
                        # store index of discareded obsCode
                        undesiredobsCodeIndex.append(k)

                    # Every 13 obsCodes is at end of line. In this case read next line and continue
                    if np.mod(k+1, 13) == 0 and nObs != 13:
                        numHeaderLines = numHeaderLines + 1
                        line = fid.readline().rstrip()
                        line = line[0:60]     # deletes 'SYS / # / OBS TYPES'
                        line_ = [el for el in line.split(" ") if el != ""]

                numOfObsCodes.append(len(GNSSSystemObsCodes[Sys]))
                obsCodes[numGNSSsystems] = GNSSSystemObsCodes
                obsCodeIndex[numGNSSsystems] = desiredObsCodeIndex # Store indices of desired obsCodes


        if 'TIME OF FIRST OBS' in line:
            line = line[0:60]     #  deletes 'TIME OF FIRST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            for k in np.arange(0,6):
                tok = line_.pop(0)
                if k ==0:
                    yyyy = int(tok)
                elif k ==1:
                    mm = int(tok)
                elif k ==2:
                    dd = int(tok)
                elif k ==3:
                    hh = int(tok)
                elif k ==4:
                    mnt = int(tok)
                elif k ==5:
                    ss = float(tok)


            tFirstObs = np.array([[yyyy],[mm],[dd],[hh],[mnt],[ss]])

            # Get Time system
            if len(line_) != 0:
                aux = line_.pop(0)
                if aux == 'GPS':
                    timeSystem = 'GPS'
                elif aux == 'GLO':
                    timeSystem = 'GLO'
                elif aux == 'GAL':
                    timeSystem = 'GAL'
                elif aux == 'BDT':
                    timeSystem = 'BDT'

            else:
                try:
                    if gnssType == 'G':
                        timeSystem = 'GPST'
                    elif gnssType == 'R':
                        timeSystem = 'GLOT'
                    elif gnssType == 'E':
                        timeSystem = 'GALT'
                    elif gnssType == 'C':
                        timeSystem = 'BDT'
                except:
                    timeSystem = "GPS"


        if 'TIME OF LAST OBS' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            for k in np.arange(0,6):
                tok = line_.pop(0)
                if k ==0:
                    yyyy = int(tok)
                elif k ==1:
                    mm = int(tok)
                elif k ==2:
                    dd = int(tok)
                elif k ==3:
                    hh = int(tok)
                elif k ==4:
                    mnt = int(tok)
                elif k ==5:
                    ss = float(tok)

            tLastObs = np.array([[yyyy],[mm],[dd],[hh],[mnt],[ss]])

        if 'INTERVAL' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            tInterval = float(line_.pop(0))


          ## -- This is an optional record!
          # if 'RCV CLOCK OFFS APPL' in line:
          #     if (strtok(line)=='0'):
          #         clockOffsetsON = 0;
          #     elif (strtok(line)=='1'):
          #         clockOffsetsON = 1;
          #     else:
          #         success = 0;
          #         print('ERROR (rinexReadsObsHeader304): unrecognized receiver clock offsets flag!')
          #         fid.close()


           ## This is an optional record
        if 'LEAP SECONDS' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            leapSec = int(line_.pop(0))


           ## -- store approximate receiver position
        if 'APPROX POSITION XYZ' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            approxPosition = np.array([[float(line_[0])],[float(line_[1])],[float(line_[2])]])


         ## GLOANSS SLOTS
        if 'GLONASS SLOT / FRQ #' in line:
            line = line[0:60]     #  deletes 'GLONASS SLOT / FRQ #'
            line_ = [el for el in line.split(" ") if el != ""]
            nGLOSat = int(line_.pop(0))
            slotNumbers = np.array([])
            channels = np.array([])
            for k in np.arange(0,nGLOSat):
                slotNumber = line_.pop(0)[1::]
                channel = int(line_.pop(0))
                slotNumbers = np.append(slotNumbers,slotNumber)
                channels = np.append(channels,channel)

                # GLONASS SLOT / FRQ # records hold up to 8 satellites per
                # line. After every 8th entry read the continuation line, but
                # only if we have not yet consumed the last satellite -- the
                # previous "k+1 == 24" check assumed exactly 24 satellites and
                # broke for partial constellations or extended ones.
                if np.mod(k + 1, 8) == 0 and (k + 1) < nGLOSat:
                    line = fid.readline().rstrip()
                    numHeaderLines = numHeaderLines + 1
                    line = line[0:60]
                    line_ = [el for el in line.split(" ") if el != ""]

            GLO_Slot2ChannelMap = dict(zip(slotNumbers.astype(int),channels.astype(int)))

        if 'REC # / TYPE / VERS' in line:
            recType = line[20:40]


     # End of Gobbling Header Loop
    for k in np.arange(0,numGNSSsystems):
        # Give info if any of GNSS systems had zero of desired obscodes.
        if numOfObsCodes[k] == 0 or sum(tFirstObs) == 0:
            if GNSSsystems[k] == 'G':
                print('INFO: (rinexReadsObsHeader304)\nNone of the GPS satellites had any of the desired obsCodes\n\n')
            elif GNSSsystems[k] == 'R':
                print('INFO: (rinexReadsObsHeader304)\nNone of the GLONASS satellites had any of the desired obsCodes\n\n')
            elif GNSSsystems[k] == 'E':
                print('INFO: (rinexReadsObsHeader304)\nNone of the Galileo satellites had any of the desired obsCodes\n\n')
            elif GNSSsystems[k] == 'C':
                print('INFO: (rinexReadsObsHeader304)\nNone of the BeiDou satellites had any of the desired obsCodes\n\n')

    ## store rinex header info
    rinexHeader['rinexVersion'] =rinexVersion
    rinexHeader['rinexType'] = rinexType
    rinexHeader['gnssType'] =gnssType
    rinexHeader['rinexProgr'] =rinexProgr
    rinexHeader['rinexDate'] =rinexDate


    print('INFO(rinexReadObsFileHeader304): Rinex header has been read')

    return success, rinexVersion, gnssType, markerName, recType, antDelta,GNSSsystems,numOfObsCodes, \
    obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval,timeSystem, numHeaderLines, clockOffsetsON, \
    rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, eof, fid


def rinexReadObsBlock304(fid, numSV, nObsCodes, GNSSsystems, obsCodeIndex, readSS, readLLI):
    """
    Reads all the observations from a RINEX observation block.

    Positioned at the beginning of the line immediately after the header of the
    observations block, reads all the observations in this block of a RINEX
    observations file. This function is meant to be used after using function
    rinexReadObsFileHeader304

    Based in the work of Antonio Pestana, rinexReadObsBlock211, March 2015
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS:
    -------

    fid:                  Matlab file identifier of a Rinex observations text file

    numSV:                total number of satellites with observations in
                          current observation block, integer

    numOfObsCodes:        column vector containing number of observation
                          types for each GNSS system. Order is the same as
                          GNSSsystems

    GNSSsystems:          cell array containing codes of GNSS systems included
                          in RINEX observationfile. Elements are strings.
                          ex. "G" or "E"

    obsCodeIndex:         cell with one cell element for each GNSS system.
                          Order is the same as GNSSsystems. Each cell element
                          contains an array of indices. These indices
                          indicate the observation types that should be
                          read for each GNSS system. ex. If one index for
                          GPS is 1 then the first observation type for GPS
                          should be read.

    readSS:                   Boolean, 0 or 1.
                              1 = read "Signal Strength" Indicators
                              0 = do not read "Signal Strength" Indicators

    readLLI:                  Boolean, 0 or 1.
                              1 = read "Loss-Of-Lock Indicators"
                              0 = do not read "Loss-Of-Lock Indicators"
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    --------

    success:               Boolean. 1 if the function seems to be successful,
                          0 otherwise

    Obs:                  matrix [numSV x max_nObs] that stores all
                          observations of this observation block. max_nObs
                          is the highest number of observation codes that
                          any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    SVlist:               column cell [numSV x 1] that conatins the
                          identification code of each line of observation
                          block. ex. "G21". numSV is total number of
                          satellites minus amount of satellites removed.

    numSV:                numSV, unlike the input of same name, is the total
                          number of satellites minus amount of satellites
                          removed.

    LLI:                  matrix [numSV x max_nObs] that stores all
                          "loss-of-lock" indicators of this observation block.
                          max_nObs is the highest number of observation codes
                          that any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    SS:                   matrix [numSV x max_nObs] that stores all
                          "signal strength" indicators of this observation block.
                          max_nObs is the highest number of observation codes
                          that any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    eof:                  end-of-file flag; 1 if end-of-file was reached,
                          0 otherwise
    --------------------------------------------------------------------------------------------------------------------------
    """
    ## Initialize variables in case of input error
    success                     = np.nan
    eof                         = np.nan
    max_n_obs_Types             = np.nan
    Obs                         = np.nan
    LLI                         = np.nan
    SS                          = np.nan
    SVlist                      = np.nan
    removed_sat                 = np.nan
    desiredGNSSsystems          = np.nan

    ## -- Testing input arguments
    if type(numSV) != int:
        print(f'INPUT ERROR(rinexReadObsBlock304): The input argument numSV is of type {type(numSV)}.\n Must be of type double')
        success = 0
        return success

    nObsCodes = [int(x) for x in nObsCodes]
    ## Test type of numOfObsCodes
    if type(nObsCodes[0]) != int:
        print('INPUT ERROR(rinexReadObsBlock304): The input argument numOfObsTypes is of type %s.\n Must be of type double' % (type(nObsCodes)))
        success = 0
        return success


    ## Test size of numOfObsCodes
    if len(nObsCodes) != len(GNSSsystems):
        print('INPUT ERROR(rinexReadObsBlock304): The input argument numOfObsTypes must have same length as GNSSsystems')
        success = 0
        return success

    success = 1
    eof     = 0

    # Highest number of obs codes of any GNSS system
    max_n_obs_Types = max(nObsCodes)

    # Initialize variables
    Obs = np.empty([numSV, max_n_obs_Types])
    SVlist = [np.nan]*numSV
    if readLLI:
        LLI = np.empty([numSV, max_n_obs_Types])

    if readSS:
        SS  = np.empty([numSV, max_n_obs_Types])

    # number of satellites excluded so far
    removed_sat = 0
    desiredGNSSsystems = list(GNSSsystems.values())
    # Gobble up observation block
    for sat in np.arange(0,numSV):
        line = fid.readline().rstrip()
        if not line:
            return

        SV = line[0:3].strip() # Satellite code, ex. 'G11' or 'E03'
        if SV[0] not in desiredGNSSsystems:
            removed_sat +=1
        else:
            ## Index of current GNSS system
            GNSSsystemIndex = [i for i in GNSSsystems if GNSSsystems[i]==SV[0]][0]
            SVlist[sat - removed_sat] = SV # Store SV of current row
            n_obs_current_system = nObsCodes[GNSSsystemIndex-1]
            for obs_num in np.arange(0, n_obs_current_system):
                obsIndex = obsCodeIndex[GNSSsystemIndex][obs_num]
                charPos = 4+(obsIndex)*16
                ## check that the current observation of the current GNSS system
                ## is not on the list of obs types to be excluded
                ## stringlength of next obs.
                obsLen = min(14, len(line) - charPos)
                # read next obs
                newObs = line[charPos:charPos+obsLen].strip()
                # If observation missing, set to 0
                if newObs != '':
                    newObs = float(newObs)
                else:
                    newObs = 0
                # Store new obs
                Obs[sat - removed_sat, obs_num] = newObs
                # read LLI of current obs (if present)
                if readLLI:
                    if charPos+13<len(line):
                        newLLI = line[charPos+13]
                    else:
                        newLLI = ' '

                    if newLLI.isspace():
                        newLLI = -999
                    else:
                        newLLI = int(newLLI)
                    LLI[sat - removed_sat, obs_num] = newLLI


                if readSS:
                    # read SS of current obs (if present)
                    if charPos+14<len(line):
                        newSS = line[charPos+14]
                    else:
                        newSS = ' '

                    # if no SS set to -999
                    if newSS.isspace():
                        newSS = -999
                    else:
                        newSS = int(newSS)
                    SS[sat - removed_sat, obs_num]  = newSS


    ## -- Update number og satellites after satellites have been excluded
    numSV = numSV - removed_sat
    ## --Remove empty arrays
    SVlist = list(filter(None,SVlist))
    idx_keep = len(Obs) -1 -removed_sat + 1 # removing sats
    Obs = Obs[:idx_keep,:]
    return success, Obs,SVlist, numSV, LLI, SS, eof


def rinexReadObsBlockHead304(fid):
    """
    Reads the metadata in the head of a RINEX 3.xx observations block, NOT
    the header of the file.

    ATTENTION: Ignores all data in blocks with event flags with numbers
    greater than 1.

    Positioned in a RINEX 3.04 GNSS observations text file at the beginning
    of an observation block. In rinex 3.xx the line starts with '> '

    --------------------------------------------------------------------------------------------------------------------------
     INPUTS

     fid:              Python identifier of an open RINEX 3.04 GNSS
                       observations text file positioned at the beginning
                       of an observation block.
    --------------------------------------------------------------------------------------------------------------------------
     OUTPUTS

     success:          1 if function performs successfully, 0 otherwise

     epochflag:        Rinex observations epoch flag, as follows:
                           0: OK
                           1: power failure between previous and current epoch
                       From now on the "event flags":
                           2: start moving antenna
                           3: new site occupation
                           4: header information follows
                           5: external event (epoch is significant)

     clockOffset:          value of the receiver clock offset. If not present
                           in the metadata of the observations block
                           (it's optional RINEX 3.04 data)it is assumed to be
                           zero. If not zero implies that epoch, code, and
                           phase data have been corrected by applying
                           realtime-derived receiver clock offset

     date:                 time stamp of the observations block. Six-elements column-vector
                           as follows:
                               year: four-digits year (eg: 1959)
                               month: integers 1..12
                               day: integers 1..31
                               hour: integers 0..24
                               minute: integers 0..60
                               second: reals 0..60

     numSV:                number of satellites with observations in with
                           observations. This will include all satellite
                           systems.
    --------------------------------------------------------------------------------------------------------------------------

    """
    # Initialize variables
    success = 1
    eof     = 0
    date    = [0,0,0,0,0,0]
    numSV   = 0
    epochflag = 0
    clockOffset = 0
    noFlag = 1

    line = fid.readline().rstrip()
    if not line:
        eof = 1
        return success, epochflag, clockOffset, date, numSV, eof
    epochflag   = line[31]
    # skip to next block if event flag is more than 1
    while int(epochflag) > 1:
        noFlag = 0
        linejump = int(line[32:35])
        msg = f'WARNING(rinexReadsObsBlockHead304): Observations event flag encountered. Flag = {str(epochflag)}, hence {str(linejump)} lines were ignored. '
        for count in np.arange(0,linejump+1):
            line = fid.readline().rstrip()
        epochflag = int(line[31]) # changed from 30 to 31 29.08.2023

    numSV = int(line[32:35])
    clockOffset = 0
    if len(line) == 56:
        clockOffset = float(line[41:56])

    # Reads the time stamp of the observations block (6 numerical values)
    date = line[1::]
    date = [float(el) for el in line[1::].split(" ") if el != ""]
    date = date[:6]

    if noFlag == 0:
        msg2 = msg + 'Epoch date: %.4d %.2d %.2d %.2d:%.2d:%6.4f' % (date[0],date[1],date[2],date[3],date[4],date[5])
        print(msg2)


    return success, epochflag, clockOffset, date, numSV, eof


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------- RINEX 2.11 ---------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------

def readRinexObs211(filename, readSS=None, readLLI=None, includeAllGNSSsystems=None,includeAllObsCodes=None, \
                    desiredGNSSsystems=None, desiredObsCodes=None, desiredObsBands=None):
    """
    Program/function to read GNSS observations in RINEX V.2 observation files
    The main core of the program is 4 functions:
                                  rinexReadObsFileHeader211
                                  rinexReadObsBlockHead211
                                  rinexReadObsBlock211
                                  rinexFindNEpochs211


    To export every parameter use this code:

    GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
         obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
         rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success = \
         readRinexObs211(rinObsFilename)

    Tip: Unwanted outputs can be ignored using _

    --------------------------------------------------------------------------------------------------------------------------
    INPUTS

    filename:                 path and name of RINEX V.2 observation file,
                              string

    readSS:                   Boolean, 0 or 1.
                              1 = read "Signal Strength" Indicators
                              0 = do not read "Signal Strength" Indicators

    readLLI:                  Boolean, 0 or 1.
                              1 = read "Loss-Of-Lock Indicators"
                              0 = do not read "Loss-Of-Lock Indicators"

    includeAllGNSSsystems:    Boolean, 0 or 1.
                              1 = include alle GNSS systems(GPS, GLONASS, Galieo, BeiDou)
                              0 = include only GNSSsystems specified in desiredGNSSsystems

    includeAllObsTypes:       Boolean, 0 or 1.
                              1 = include all valid ObsTypes
                              0 = include only ObsTypes specified in desiredObsTypes

    desiredGNSSsystems:       array og strings containing  codes of desired
                              GNSSsystems to be included,
                              ex. ["G", "E"]
                              OBS: Must be string array, NOT char vector

    desiredObsTypes:          array of strings containing desired ObsTypes to be
                              included, ex. ["C", "L", "S", "D"]
                              OBS: Must be string array, NOT char vector

    desiredObsBands:          array of desired obs Bands to be included,
                              ex [1, 5]
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS

    GNSS_obs:                 cell containing a matrix for each GNSS system.
                              Each matrix is a 3D matrix containing all
                              observation of current GNSS system for all epochs.
                              Order of obsType index is same order as in
                              obsCodes cell

                              GNSS_obs{GNSSsystemIndex}(PRN, obsType, epoch)
                                              GNSSsystemIndex: double,
                                              1,2,...,numGNSSsystems
                                              PRN: double
                                              ObsType: double: 1,2,...,numObsTypes
                                              epoch: double

    GNSS_LLI:                 cell containing a matrix for each GNSS system
                              Each matrix stores loss of lock indicators for
                              each epoch for that GNSS system

    GNSS_SS:                  cell containing a matrix for each GNSS system
                              Each matrix stores signal strength indicators
                              for each epoch for that GNSS system

    GNSS_SVs:                 cell containing a matrix for each GNSS system.
                              Each matrix contains number of satellites with
                              obsevations for each epoch, and PRN for those
                              satellites

                              GNSS_SVs{GNSSsystemIndex}(epoch, j)
                                              j=1: number of observed satellites
                                              j>1: PRN of observed satellites

    time_epochs:              matrix conatining gps-week and "time of week"
                              for each epoch
                              time_epochs(epoch,i),   i=1: week
                                                      i=2: time-of-week in seconds (tow)

    nepochs:                  number of epochs with observations in rinex observation file.

    GNSSsystems:              cell array containing codes of GNSS systems included
                              in RINEX observationfile. Elements are strings.
                              ex. "G" or "E"

    obsCodes:                 Cell that defines the observation
                              codes available for all GNSS system. Each cell
                              element is another cell containing the codes for
                              that GNSS system. Each element in this cell is a
                              string with three-characters. The first
                              character (a capital letter) is an observation code
                              ex. "L" or "C". The second character (a digit)
                              is a frequency code. The third character(a Capital letter)
                              is the attribute, ex. "P" or "X"

    approxPosition:           array containing approximate position from rinex
                              observation file header. [X, Y, Z]

    max_sat:                  array conataining max PRN number for each GNSS
                              system. Follows same order as GNSSsystems

    tInterval:                observations interval; seconds.

    markerName:               name of the antenna marker; '' if not specified

    rinexVersion:             string. rinex observation file version

    recType:                  receiver type, char vector

    timeSystem:               three-character code string of the time system
                              used for expressing tfirstObs;
                              can be GPS, GLO or GAL;

    leapSec:                  number of leap seconds since 6-Jan-1980.
                              UTC=GPST-leapSec. NaN by default.
                              THIS IS RINEX V.2 OPTIONAL DATA

                              rinexHeader: cell column-vector containing the
                              following data:
                              rinexVersion:   RINEX version number; string.
                                              '' if not specified
                              rinexType:      RINEX file type; char

    gnssType:                 GNSS system of the satellites observed; can be
                              'G', 'R', 'E', 'C' or 'M' that stand for
                              GPS, GLONASS, GALILEO, BeiDou or Mixed ; char

    rinexProgr:               name of the software used to produce de RINEX
                              GPS obs file; '' if not specified


    rinexDate:                date/time of the RINEX file creation; '' if not
                              specified; char


    antDelta:                 column vector ot the three components of the
                              distance from the marker to the antenna,
                              in the following order - up, east and north;
                              null vector by default

    tFirstObs:                time stamp of the first observation record in the RINEX
                              observations file; column vector
                              [YYYY; MM; DD; hh; mm; ss.sssssss];
                              THIS IS CRITICAL DATA

    tLastObs:                 time stamp of the last observation record in the RINEX
                              observations file; column vector
                              [YYYY, MM, DD, hh, mm,ss.sssssss]. NaN by default.
                              THIS IS RINEX V.2 OPTIONAL DATA

    clockOffsetsON:           receiver clock offsets flag. O if no realtime-derived
                              receiver clock offset was applied to epoch,
                              code and phase data (in other words, if the
                              file only has raw data), 1 otherwise. 0 by default.
                              THIS IS RINEX V.2 OPTIONAL DATA

    GLO_Slot2ChannelMap:      map container that maps GLONASS slot numbers to
                              their respective channel number.
                              GLO_Slot2ChannelMap(slotnumber)

    success:                  Boolean. 1 if the reading of the RINEX
                              observations file seems to be successful,
                              0 otherwise
    --------------------------------------------------------------------------------------------------------------------------

    ADVICE: The function rinexFindNEpochs() calculates the amount of observation epochs in
    advance. It is recommended to manually add this information to the header if it is
    not included by default.
    --------------------------------------------------------------------------------------------------------------------------

    According to RINEX V.2 the observation type codes are:

    Observation code
            C: Pseudorange
    GPS: C/A, L2C
    Glonass: C/A
    Galileo: All
    L: Carrier phase
    D: Doppler frequency
    S: Raw signal strengths or SNR values as given by the receiver for the
    respective phase observations
    I: Ionosphere phase delay
    X: Receiver channel numbers

    Frequency code
    GPS Glonass Galileo SBAS
    1: L1 G1 E1 B1    (GPS,QZSS,SBAS,BDS)
    2: L2 G2 B1-2     (GLONASS)
    4: G1a            (Galileo)
    5: L5 E5a B2/B2a  (GPS, QZSS, SBAS, IRNSS)
    6: L6 E6 B3 G2a   (Galileo, QZSS, BDS, GLONASS)
    7: E5b B2/B2b     (Galileo)
    8: E5a+b E5a+b    (Galileo, BDS)
    9: S              (IRNSS)
    0: for type X     (all)

    Attribute:
    A = A channel     (Galileo,IRNSS,GLONASS)
    B = B channel     (Galileo,IRNSS,GLONASS)
    C = C channel     (Galiloe, IRNSS)
    C code-based  (SBAS, GPS, GLONASS, QZSS)
    D = Semi-codelss  (GPS)

    I = I channel     (GPS, Galileo, QZSS, BDS)
    L = L channel     (L2C GPS, QZSS)
    P channel     (GPS. QZSS)
    M = M code-based  (GPS)
    N = Codeless      (GPS)
    P = P code-based  (GPS, GLONASS)
    Pilot channel (BDS)

    Q = Q channel     (GPS, Galileo, QZSS, BDS)
    S = D channel     (GPS, Galileo, QZSS, BDS)
    M channel     (L2C GPS, QZSS)

    W = Based on Z-tracking (GPS)
    X = B+C channels  (Galileo, IRNSS)
    I+Q channels  (GPS, IRNSS)
    M+L channels  (GPS, QZSS)
    D+P channels  (GPS, QZSS, BDS)

    Y = Y code based  (GPS)
    Z = A+B+C channels(Galileo)
    D+P channels  (BDS)
    --------------------------------------------------------------------------------------------------------------------------
    """
    ## -- Setting None arguments
    if readSS is None:
        readSS = 1
    if readLLI is None:
        readLLI = 1
    if includeAllGNSSsystems is None and desiredGNSSsystems is None:
        includeAllGNSSsystems = 1
    if includeAllObsCodes is None and desiredObsCodes is None:
        includeAllObsCodes = 1
    if desiredGNSSsystems is None:
        desiredGNSSsystems = ['G','R','E','C']
    if desiredObsCodes is None:
        desiredObsCodes = ['C','L','S','D']
    if desiredObsBands is None:
        desiredObsBands = list(np.arange(1,10))

    ## Get the start time
    t = time.process_time()

    ## - Initialize variables in case of input error
    GNSS_obs       = np.nan
    GNSS_LLI       = np.nan
    GNSS_SS        = np.nan
    GNSS_SVs       = np.nan
    time_epochs    = np.nan
    nepochs        = np.nan
    GNSSsystems    = np.nan
    obsCodes       = np.nan
    approxPosition = np.nan
    max_sat        = np.nan
    tInterval      = np.nan
    markerName     = np.nan
    rinexVersion   = np.nan
    recType        = np.nan
    timeSystem     = np.nan
    leapSec        = np.nan
    gnssType       = np.nan
    rinexProgr     = np.nan
    rinexDate      = np.nan
    antDelta       = np.nan
    tFirstObs      = np.nan
    tLastObs       = np.nan
    clockOffsetsON = np.nan
    GLO_Slot2ChannelMap = np.nan

    ### --- Dict for storing data
    GNSS_obs = {}
    GPS = {}
    GLONASS = {}
    Galileo = {}
    BeiDou = {}
    GPS_LLI = {}
    GLONASS_LLI = {}
    Galileo_LLI = {}
    BeiDou_LLI = {}
    GPS_SS = {}
    GLONASS_SS = {}
    Galileo_SS = {}
    BeiDou_SS = {}

    ## -- Test if readSS is boolean
    if readSS!=1 and readSS!=0:
        print('INPUT ERROR(readRinexObs211): The input argument readSS must be either 1 or 0')
        success = 0
        return


    ## -- Test if readLLI is boolean
    if readLLI!=1 and readLLI!=0:
        print('INPUT ERROR(readRinexObs211): The input argument readLLI must be either 1 or 0')
        success = 0
        return

    max_GPS_PRN     = 36 # Max number of GPS PRN in constellation
    max_GLONASS_PRN = 36 # Max number of GLONASS PRN in constellation
    max_Galileo_PRN = 36 # Max number of Galileo PRN in constellation
    max_Beidou_PRN  = 100 # Max number of BeiDou PRN in constellation

    #  Read header of observation file
    [success, rinexVersion, gnssType, markerName, recType, antDelta,\
    GNSSsystems,numOfObsCodes, obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval, \
    timeSystem, _, clockOffsetsON, rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, _, fid] = \
    rinexReadObsFileHeader211(filename, includeAllGNSSsystems, includeAllObsCodes,desiredGNSSsystems, desiredObsCodes, desiredObsBands)

    if success==0:
        return

    ## -- Read all remaining observation data at once for fast processing
    remaining_lines = fid.readlines()
    fid.close()

    ## -- Count epochs from buffered lines
    _epoch_re = re.compile(r"\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}\.\d+)\s*")
    epoch_line_indices = [i for i, ln in enumerate(remaining_lines) if _epoch_re.match(ln)]
    nepochs = len(epoch_line_indices)

    if nepochs == 0:
        print('ERROR(readRinexObs211): No epoch records found in observation file')
        success = 0
        return

    # Compute tInterval from first two epoch headers if not in header
    if np.isnan(tInterval) and nepochs >= 2:
        m1 = _epoch_re.match(remaining_lines[epoch_line_indices[0]])
        m2 = _epoch_re.match(remaining_lines[epoch_line_indices[1]])
        if m1 and m2:
            tInterval = float(m2.group(6)) - float(m1.group(6))

    # Compute tLastObs if not in header
    try:
        _tl_is_nan = np.all(np.isnan(tLastObs))
    except (TypeError, ValueError):
        _tl_is_nan = False
    if _tl_is_nan:
        m_last = _epoch_re.match(remaining_lines[epoch_line_indices[-1]])
        if m_last:
            yr = float(m_last.group(1))
            try:
                full_year = float(str(int(tFirstObs[0][0]))[0:2] + str(int(yr)))
            except Exception:
                full_year = float('20' + str(int(yr)))
            tLastObs = np.array([[full_year], [float(m_last.group(2))], [float(m_last.group(3))],
                                [float(m_last.group(4))], [float(m_last.group(5))], [float(m_last.group(6))]])
        print('INFO(readRinexObs211): TIME OF LAST OBS computed from observation data')

    success = 1

    ## -- Create a file-like cursor over pre-read lines for block parsers
    class _LineCursor:
        __slots__ = ('_lines', '_pos', '_n')
        def __init__(self, lines):
            self._lines = lines
            self._pos = 0
            self._n = len(lines)
        def readline(self):
            if self._pos >= self._n:
                return ''
            line = self._lines[self._pos]
            self._pos += 1
            return line
    _cursor = _LineCursor(remaining_lines)

    # Number of GNSS systems
    nGNSSsystems = len(GNSSsystems)

    # Declare data cells, arrays and matrices
    GNSS_SVs = {}
    max_sat  =  np.zeros([nGNSSsystems,1])
    t_week = []
    t_tow = []

    GNSSsystems_full_names =  [""]*nGNSSsystems
    GNSS_LLI = {}
    GNSS_SS = {}

    # Create array for max_sat. Initialize cell elements in cell arrays
    for k in range(nGNSSsystems):
        if GNSSsystems[k + 1] == 'G':
            max_sat[k] = max_GPS_PRN
            GNSS_SVs['G'] = np.zeros([nepochs, max_GPS_PRN + 1], dtype=np.int16)
            GNSSsystems_full_names[k] = "GPS"
        elif GNSSsystems[k + 1] == 'R':
            max_sat[k] = max_GLONASS_PRN
            GNSS_SVs['R'] = np.zeros([nepochs, max_GLONASS_PRN + 1], dtype=np.int16)
            GNSSsystems_full_names[k] = "GLONASS"
        elif GNSSsystems[k + 1] == 'E':
            max_sat[k] = max_Galileo_PRN
            GNSS_SVs['E'] = np.zeros([nepochs, max_Galileo_PRN + 1], dtype=np.int16)
            GNSSsystems_full_names[k] = "Galileo"
        elif GNSSsystems[k + 1] == 'C':
            max_sat[k] = max_Beidou_PRN
            GNSS_SVs['C'] = np.zeros([nepochs, max_Beidou_PRN + 1], dtype=np.int16)
            GNSSsystems_full_names[k] = "BeiDou"
        else:
            print(f'ERROR(readRinexObs211): Only following GNSS systems are compatible with this program: GPS, GLONASS, Galileo, Beidou. {GNSSsystems[k]} is not valid')
            return

        curr_sys = GNSSsystems[k+1]
        # Note: GNSS_obs[curr_sys] is populated as a per-epoch dict below
        # (e.g. GNSS_obs['G'] = GPS), so no large 3D pre-allocation is needed here.

        # Preallocation LLI and SS
        if readLLI:
            GNSS_LLI[curr_sys] = np.zeros([nGNSSsystems,1])
        else:
            GNSS_LLI[curr_sys] = np.nan
        if readSS:
            GNSS_SS[curr_sys] = np.zeros([nGNSSsystems,1])
        else:
            GNSS_SS[curr_sys] = np.nan

    GNSS_names = dict(zip(['G', 'R', 'E', 'C'],['GPS','GLONASS','Galileo','Beidou']))

    ## -- Pre-compute system-to-index lookup
    sys_to_idx = {v: k for k, v in GNSSsystems.items()}
    max_sat_int = {k+1: int(max_sat[k].item()) if hasattr(max_sat[k], "item") else int(max_sat[k]) for k in range(nGNSSsystems)}

    ## -- Per-system storage dicts
    _obs_dicts = {'G': GPS, 'R': GLONASS, 'E': Galileo, 'C': BeiDou}
    _lli_dicts = {'G': GPS_LLI, 'R': GLONASS_LLI, 'E': Galileo_LLI, 'C': BeiDou_LLI}
    _ss_dicts  = {'G': GPS_SS, 'R': GLONASS_SS, 'E': Galileo_SS, 'C': BeiDou_SS}

    current_epoch = 0
    n_update_break = max(int(nepochs // 10), 1)
    bar_format = '{desc}: {percentage:3.0f}%|{bar}| ({n_fmt}/{total_fmt})'

    with tqdm(total=100, desc="Rinex observations are being read", position=0, leave=True, bar_format=bar_format) as pbar:
        while 1:
            #  Read Obs Block Header
            success, _, _, date, numSV, SVlist_, eof = rinexReadObsBlockHead211(_cursor)
            if success==0 or eof==1:
                break
            # Read current block of observations
            success, Obs, SVlist, numSV, LLI, SS, eof = rinexReadObsBlock211(_cursor, numSV, numOfObsCodes, GNSSsystems, obsCodeIndex, readSS, readLLI, SVlist_)
            if success==0 or eof==1:
                break

            current_epoch += 1
            if current_epoch % n_update_break == 0:
                pbar.update(10)

            ## Convert date to GPS-week and "time-of-week"
            date[0] = float(str(tFirstObs[0][0])[0:2] + str(int(date[0])))
            week, tow = date2gpstime(int(date[0]), int(date[1]), int(date[2]), int(date[3]), int(date[4]), int(date[5]))
            t_week.append(week)
            t_tow.append(tow)

            ## Initialize per-epoch arrays
            nGNSS_sat_current_epoch = [0] * nGNSSsystems
            GNSS_obs_dum = {}
            GNSS_LLI_dum = {}
            GNSS_SS_dum  = {}
            for k in range(nGNSSsystems):
                mk = max_sat_int[k+1]
                nc = numOfObsCodes[k]
                GNSS_obs_dum[k+1] = np.zeros((mk + 1, nc))
                GNSS_LLI_dum[k+1] = np.zeros((mk + 1, nc))
                GNSS_SS_dum[k+1]  = np.zeros((mk + 1, nc))

            ## -- Iterate through satellites of epoch and store obs, LLI and SS
            for sat in range(numSV):
                curr_sys = SVlist[sat][0]
                gi = sys_to_idx[curr_sys]
                nGNSS_sat_current_epoch[gi - 1] += 1
                SV = int(SVlist[sat][1:3])
                n_obs = numOfObsCodes[gi - 1]
                GNSS_obs_dum[gi][SV][:n_obs] = Obs[sat, :n_obs]
                if readLLI:
                    GNSS_LLI_dum[gi][SV][:n_obs] = LLI[sat, :n_obs]
                if readSS:
                    GNSS_SS_dum[gi][SV][:n_obs] = SS[sat, :n_obs]
                GNSS_SVs[curr_sys][current_epoch-1, nGNSS_sat_current_epoch[gi-1]] = SV

            for k in range(nGNSSsystems):
                curr_sys = GNSSsystems[k+1]
                GNSS_SVs[curr_sys][current_epoch-1, 0] = nGNSS_sat_current_epoch[k]
                _obs_dicts[curr_sys][current_epoch] = GNSS_obs_dum[k+1]
                if readLLI:
                    _lli_dicts[curr_sys][current_epoch] = GNSS_LLI_dum[k+1]
                if readSS:
                    _ss_dicts[curr_sys][current_epoch] = GNSS_SS_dum[k+1]

        ## -- Build time_epochs array (once, at end)
        time_epochs = np.column_stack((t_week, t_tow)) if t_week else np.empty((0, 2))

        ## -- Storing observation in dictionary
        GNSS_obs['G'] = GPS
        GNSS_obs['R'] = GLONASS
        GNSS_obs['E'] = Galileo
        GNSS_obs['C'] = BeiDou
        ## -- Storing "loss of lock indicaors"  in dict
        GNSS_LLI['G'] = GPS_LLI
        GNSS_LLI['R'] = GLONASS_LLI
        GNSS_LLI['E'] = Galileo_LLI
        GNSS_LLI['C'] = BeiDou_LLI
        ## -- Storing  SS in dict
        GNSS_SS['G'] = GPS_SS
        GNSS_SS['R'] = GLONASS_SS
        GNSS_SS['E'] = Galileo_SS
        GNSS_SS['C'] = BeiDou_SS

        del_sys = list(GNSS_obs.keys())
        for sys in del_sys: # Deleting systems with no observations
            if not GNSS_obs[sys]:
                del GNSS_obs[sys]

        ## -- Removing system not in desiredGNSSsystems
        del_GNSS_LLI = [k for k,v in GNSS_LLI.items() if k not in desiredGNSSsystems]
        del_GNSS_SS  = [k for k,v in GNSS_SS.items() if k not in desiredGNSSsystems]
        del_GNSS_obs = [k for k,v in GNSS_obs.items() if k not in desiredGNSSsystems]
        del_GNSS_names = [k for k,v in GNSS_names.items() if k not in desiredGNSSsystems]
        del_GNSSsystems = [k for k,v in GNSSsystems.items() if v not in desiredGNSSsystems]
        del_GNSS_SVs = [k for k,v in GNSS_SVs.items() if k not in desiredGNSSsystems]

        keep_idx = []
        for k, v in GNSSsystems.items():
            if v in desiredGNSSsystems:
                keep_idx.append(k-1)

        max_sat = max_sat[keep_idx]

        for k in del_GNSS_LLI:
            del GNSS_LLI[k]
        for k in del_GNSS_SS:
            del GNSS_SS[k]
        for k in del_GNSS_obs:
            del GNSS_obs[k]
        for k in del_GNSS_names:
            del GNSS_names[k]
        for k in del_GNSSsystems:
            del GNSSsystems[k]
        for k in del_GNSS_SVs:
            del GNSS_SVs[k]


        GNSSsystems = {i+1: v for i, v in enumerate(GNSSsystems.values())} #update keys

        if current_epoch!= nepochs and success == 1:
            print('ERROR(readRinexObs211): The amount of epochs calculated in advance(nepochs = %d) does not equal number og epochs prossesed(current_epoch = %d).\nCheck that header information concerning TIME OF FIRST OBS and TIME OF LAST OBS is correct.\n' %(nepochs, current_epoch))


        messages = {}
        if success == 1:
            messages[0]= 'INFO(readRinexObs211): The following GNSS systems have been read into the data:'
            for curr_sys in list(GNSSsystems.values()):
                k = list(GNSSsystems.values()).index(curr_sys)
                messages[k+1]= 'INFO(readRinexObs211): The following %s observation types have been registered:' % (GNSS_names[curr_sys])
                for obs in np.arange(0, len(obsCodes[k+1][curr_sys])):
                    if obs == 0:
                        messages[k+1]= messages[k+1] + ' %s' % (obsCodes[k+1][curr_sys][obs])
                    else:
                        messages[k+1]= messages[k+1] + ', %s' % (obsCodes[k+1][curr_sys][obs])

                if k == 0:
                    messages[0]= messages[0] + ' %s' % GNSS_names[GNSSsystems[1]]
                else:
                    messages[0]= messages[0] + ', %s' % GNSS_names[GNSSsystems[k+1]]

            for msg in np.arange(0,len(messages)):
                print(messages[msg])


        if readLLI:
            print('INFO(readRinexObs211): LLI have been read (if present in observation file)')
        else:
            print('INFO(readRinexObs211): LLI have not been read')


        if readSS:
            print('INFO(readRinexObs211): SS have been read (if present in observation file)')
        else:
            print('INFO(readRinexObs211): SS have not been read')


        ## --  Finding processing time
        et = time.process_time()  # get the end time
        e = et - t                # get execution time

        if e >= 3600:
            hours = np.floor(e/3600)
            minutes = np.floor((e-hours*3600)/60)
            seconds = e-hours*3600-minutes*60
            print('INFO(readRinexObs211): Total processing time: %d hours, %d minutes, %f seconds\n' % (hours, minutes, seconds))
        elif e>60:
            minutes = np.floor(e/60)
            seconds = e-minutes*60
            print('INFO(readRinexObs211): Total processing time: %d minutes, %f seconds\n' % (minutes, seconds))
        else:
            print('INFO(readRinexObs211): Total processing time: %f seconds\n\n' % (e))


    return GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems,\
        obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType,\
        rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success


def rinexFindNEpochs211(filename, tFirstObs, tLastObs, tInterval):
    """
    Function that computes number of epochs in Rinex 3.xx observation file.
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS

    filename:         RINEX observation filename

    tFirstObs:        time stamp of the first observation record in the RINEX
                      observations file; column vector
                      [YYYY; MM; DD; hh; mm; ss.sssssss];

    tLastObs:         time stamp of the last observation record in the RINEX
                      observations file; column vector
                      [YYYY; MM; DD; hh; mm; ss.sssssss]. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function

    tInterval:        observations interval; seconds. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function.
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    -------

    nepochs:          number of epochs in Rinex observation file with
                      observations

    tLastObs:         time stamp of the last observation record in the RINEX
                      observations file; column vector
                      [YYYY, MM, DD, hh, mm, ss.sssssss]. If this information
                      was not available in rinex observation header the
                      default value is Nan. In this case the variable is
                      determined in this function

    tInterval:        observations interval; seconds. If this information
                      was not available in rinex observation header the
                      default value is Nan.

    success:                  Boolean. 1 if the function seems to be successful,
                              0 otherwise
    --------------------------------------------------------------------------------------------------------------------------

    ADVICE: The function rinexFindNEpochs() calculates the amount of observation epochs in
    advance. This calculation will be incredibly more effective if TIME OF
    LAST OBS is included in the header of the observation file. It is
    strongly advized to manually add this information to the header if it is
    not included by default.
    --------------------------------------------------------------------------------------------------------------------------
    """
    # #  Testing input arguments
    success = 1
    nepochs = 0

    ## --Test if filename is valid format
    if type(filename) is not str:
        raise TypeError('INPUT ERROR(rinexFindNEpoch): The input argument filename'\
            'is of type %s. Must be of type string' % type(filename))


    ## --Open observation file
    fid = open(filename, 'rt')
    seconds_in_a_week = 604800
    #  tLastObs is in header
    if ~np.all(np.isnan(tLastObs)):  # endret 07.12.2022
        #  tInterval is not in header
        if np.isnan(tInterval):
            tInterval_found = 0
            first_epoch_found = 0
            #  calculate tInterval
            while not tInterval_found: #  calculate tInterval
                line = fid.readline().rstrip()
                #  start of new epoch
                pattern = r"\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}\.\d+)\s*"
                if re.match(pattern, line):
                # if '>' in line:
                    if not first_epoch_found: #  first epoch
                        first_epoch_time = [el for el in line.split(" ") if el != ""][0:6]
                        first_epoch_time = [float(el) for el in first_epoch_time]
                        first_epoch_found = 1
                    else: #  seconds epoch
                        second_epoch_time = [el for el in line.split(" ") if el != ""][0:6]
                        second_epoch_time = [float(el) for el in second_epoch_time]
                        tInterval = second_epoch_time[5]-first_epoch_time[5]
                        tInterval_found = 1

        # fid.close(); fid = open(filename, 'rt')
        file = open(filename, 'rt')
        tFirstObs = tFirstObs.astype(int)
        tLastObs = tLastObs.astype(int)
        rinex_lines = file.readlines()
        # idx_start = rinex_lines.index('                                                            END OF HEADER       \n')
        idx_start = [i for i, line in enumerate(rinex_lines) if line.strip() == "END OF HEADER"][0]
        rinex_lines = rinex_lines[idx_start::]
        pattern = r"\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}\.\d+)\s*"
        epoch_line = [line for line in rinex_lines if re.search(pattern, line)] # list with all the line thats defines a epoch
        nepochs = len(epoch_line)
        file.close()
    #  if tLastObs is not in header. Function counts number of epochs manually
    else:
    ## New code for finding last
        print('INFO(rinexFindEpochs211): The header of the rinex observation file does not contain TIME OF LAST OBS.\n' \
            'This will be calculated, but consider editing rinex header to include TIME OF LAST HEADER')

        ## -- Find tLastObs
        pattern = r"\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}\.\d+)\s*"
        tLastObs = find_match_in_file(filename, pattern)
        tLastObs = np.array([list(map(float, tLastObs))]).T
        try:
            tLastObs[0] = tFirstObs[0]
        except:
            tLastObs[0] = float('20' + str(tLastObs[0][0]))
            print("\n Neither tFirstObs or tLastObs exist in header. The year is therefor a guess for after\
                  year 2000")

        first_ep,second_ep = find_first_two_epochs(filename, pattern)

        ## Computing the tInterval if not present in the header
        if np.isnan(tInterval):
            first_epoch_time  = [float(el) for el in first_ep]
            second_epoch_time = [float(el) for el in second_ep]
            tInterval = second_epoch_time[5]-first_epoch_time[5]

        # nepochs = time_difference(tFirstObs, tLastObs)/tInterval
        # fid.close(); fid = open(filename, 'rt')
        file = open(filename, 'rt')
        rinex_lines = file.readlines()
        idx_start = [i for i, line in enumerate(rinex_lines) if line.strip() == "END OF HEADER"][0]
        rinex_lines = rinex_lines[idx_start::]
        pattern = r"\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}\.\d+)\s*"
        epoch_line = [line for line in rinex_lines if re.search(pattern, line)] # list with all the line thats defines a epoch
        nepochs = len(epoch_line)
        print('INFO(rinexFindNEpochs211): TIME OF LAST OBS has been found and amount of epochs have been computed')
        file.close()
    return int(nepochs), tLastObs, tInterval, success


def find_match_in_file(file_path, pattern):
    """
    Function that reads in a file backwards and looks for pattern from the
    bottom and up. Is used for finding tLastObs when not defined in header.
    Uses binary mode to save processing time.
    """
    with open(file_path, "rb") as f:
        # Go to the end of the file
        f.seek(0, 2)
        file_size = f.tell()
        # Read the file backwards
        while file_size > 0:
            # Read a chunk of the file
            chunk_size = min(1024, file_size)
            f.seek(file_size - chunk_size)
            chunk = f.read(chunk_size)
            # Split the chunk into lines
            lines = chunk.split(b"\n")
            # Process the lines in reverse order
            for line in reversed(lines):
                line = line.decode("utf-8").rstrip()
                match = re.match(pattern, line)
                if match:
                    return match.groups()
            # Repeat until the entire file has been processed
            file_size -= chunk_size
    # No match was found
    return None


def find_match_in_file_FWD(file_path, pattern):
    """
    Function that reads in a file and looks for pattern from the bottom and up.
    Is used for finding tLastObs when not defined in header. Uses binary mode
    to save processing time. Returns current line and next 5 lines after match.
    """
    with open(file_path, "rb") as f:
        # Go to the end of the file
        f.seek(0, 2)
        file_size = f.tell()
        # Read the file backwards
        while file_size > 0:
            # Read a chunk of the file
            chunk_size = min(1024, file_size)
            f.seek(file_size - chunk_size)
            chunk = f.read(chunk_size)
            # Split the chunk into lines
            lines = chunk.split(b"\n")
            # Process the lines in reverse order
            for i in range(len(lines)-1, -1, -1):
                line = lines[i].decode("utf-8").rstrip()
                match = re.match(pattern, line)
                if match:
                    return [line] + [lines[j].decode("utf-8").rstrip() for j in range(i+1, min(i+6, len(lines)))]
            # Repeat until the entire file has been processed
            file_size -= chunk_size
    # No match was found
    return None


def find_first_two_epochs(file_path, pattern):
    """
    Function that extracts the two first epochs for
    RINEX obs file to compute tInterval when not defined.
    """
    with open(file_path, 'r') as file:
        found_header = False
        count = 0
        matches = []
        for line in file:
            if not found_header:
                if line.strip() == "END OF HEADER":
                    found_header = True
                continue
            match = re.search(pattern, line)
            if match:
                count += 1
                matches.append(match.groups())
                if count >= 2:
                    break
        ep1, ep2  = matches
        ep1 = list(ep1)
        ep2 = list(ep2)
        return ep1,ep2


def find_nepochs(file_path, pattern):
    """Find number of epochs"""
    with open(file_path, 'r') as file:
        contents = file.read()
        return len(re.findall(pattern, contents))


def time_difference(arr1, arr2):
    """Find time difference between two arrays"""
    time1 = datetime(year=int(arr1[0]), month=int(arr1[1]), day=int(arr1[2]), hour=int(arr1[3]), minute=int(arr1[4]), second=int(arr1[5]))
    time2 = datetime(year=int(arr2[0]), month=int(arr2[1]), day=int(arr2[2]), hour=int(arr2[3]), minute=int(arr2[4]), second=int(arr2[5]))
    time_delta = (time2-time1).total_seconds()
    return time_delta

def rinexReadObsFileHeader211(filename, includeAllGNSSsystems, includeAllObsCodes,desiredGNSSsystems,desiredObsCodes, desiredObsBands):
    """
    Extracts relevant data from the header of a RINEX 3.xx GNSS observations
    file. Excludes undesired GNSS systems, obsevation codes and/or frequency
    bands.

    --------------------------------------------------------------------------------------------------------------------------
    INPUTS:
    ------

    filename:                     RINEX observation filename and path

    includeAllGNSSsystems:        Boolean, 0 or 1.
                                      1 = include alle GNSS systems
                                          (GPS, GLONASS, Galieo, BeiDou)
                                      0 = include only GNSSsystems
                                          specified in desiredGNSSsystems

    includeAllobsCodes:           Boolean, 0 or 1.
                                      1 = include all valid obsCodes
                                      0 = include only obsCodes
                                          specified in desiredobsCodes

    desiredGNSSsystems:           string array containing desired GNSSsystems
                                  to be included, ex. ["G", "E", "C"]

    desiredobsCodes:              string array containing desired obsCodes to
                                  be included, ex. ["C", "L", "S", "D"]

    desiredObsBands:              array of desired obs Bands to be included,
                                  ex [1, 5]

    NOTE: If both includeAllGNSSsystems and includeAllobsCodes Boolean are 1
          then the last three input arguments are optional to include and may
          be left blank without en error.
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    -------

    success:                      1 if the reading of the RINEX observations
                                  file seems to be successful, 0 otherwise

    rinexVersion:                 string. rinex observation file version

    gnssType:                     GNSS system of the satellites observed; can
                                  be 'G', 'R', 'E', 'C' or 'M' that stand for
                                  GPS, GLONASS, GALILEO, BeiDou or Mixed; char

    markerName:                   name of the antenna marker; '' if not
                                  specified

    recType:                      Receiver type, char vector

    antDelta:                     column vector ot the three components of
                                  the distance from the marker to the antenna,
                                  in the following order - up, east and north;
                                  null vector by default

    GNSSsystems:                  cell array containing codes of GNSS systems
                                  included in RINEX observationfile. Elements
                                  are strings. ex. "G" or "E"

    numOfObsCodes:                column vector containing number of observation
                                  types for each GNSS system. Order is the same
                                  as GNSSsystems

    obsCodes:                     Cell that defines the observation
                                  codes available for all GNSS system. Each
                                  cell element is another cell containing the
                                  codes for that GNSS system. Each element in
                                  this cell is a string with three-characters.
                                  The first character (a capital letter) is
                                  an observation code ex. "L" or "C". The
                                  second character (a digit) is a frequency
                                  code. The third character(a Capital letter)
                                  is the attribute, ex. "P" or "X"

    obsCodeIndex:                 cell with one cell element for each GNSS
                                  system. Order is the same as GNSSsystems.
                                  Each cell element contains an array of
                                  indices. These indices indicate the
                                  observation types that should be read
                                  for each GNSS system. ex. If one index for
                                  GPS is 1 then the first observation type
                                  for GPS should  be read.

    tFirstObs:                    time stamp of the first observation record
                                  in the RINEX observations file; column vector
                                  [YYYY; MM; DD; hh; mm; ss.sssssss];
                                  THIS IS CRITICAL DATA

    tLastObs:                     time stamp of the last observation record
                                  in the RINEX observations file; column vector
                                  [YYYY; MM; DD; hh; mm;ss.sssssss].
                                  NaN by default.
                                  THIS IS RINEX V.2 OPTIONAL DATA

    tInterval:                    observations interval; seconds.

    timeSystem:                   three-character code string of the time
                                  system used for expressing tfirstObs;
                                  can be GPS, GLO or GAL;

    numHeaderLines:               number of lines in header

    rinexProgr:                   name of the software used to produce de
                                  RINEX GPS obs file; '' if not specified

    rinexDate:                    date/time of the RINEX file creation; ''
                                  if not specified; char

    leapSec:                      number of leap seconds since 6-Jan-1980.
                                  UTC=GPST-leapSec. NaN by default.
                                  THIS IS RINEX V.2 OPTIONAL DATA

    approxPosition:               array containing approximate position from
                                  rinex observation file header. [X, Y, Z]

    GLO_Slot2ChannelMap:          map container that maps GLONASS slot
                                  numbers to their respective channel number.
                                  GLO_Slot2ChannelMap(slotnumber)

    eof:                          end-of-file flag; 1 if end-of-file was reached,
                                  0 otherwise

    fid:                          Matlab file identifier of a Rinex
                                  observations text file
    --------------------------------------------------------------------------------------------------------------------------

    According to RINEX V.2 these codes are:

       Observation codes:
       ------------------
       C: Pseudorange
          GPS: C/A, L2C
          Glonass: C/A
          Galileo: All
       L: Carrier phase
       D: Doppler frequency
       S: Raw signal strengths or SNR values as given by the receiver for the
          respective phase observations
       I: Ionosphere phase delay
       X: Receiver channel numbers

       Frequency code:
       ---------------
       GPS Glonass Galileo SBAS
       1: L1 G1 E1 B1    (GPS,QZSS,SBAS,BDS)
       2: L2 G2 B1-2     (GLONASS)
       4: G1a            (Galileo)
       5: L5 E5a B2/B2a  (GPS, QZSS, SBAS, IRNSS)
       6: L6 E6 B3 G2a   (Galileo, QZSS, BDS, GLONASS)
       7: E5b B2/B2b     (Galileo)
       8: E5a+b E5a+b    (Galileo, BDS)
       9: S              (IRNSS)
       0: for type X     (all)

       Attribute:
       ----------
       A = A channel     (Galileo,IRNSS,GLONASS)
       B = B channel     (Galileo,IRNSS,GLONASS)
       C = C channel     (Galiloe, IRNSS)
           C code-based  (SBAS, GPS, GLONASS, QZSS)
       D = Semi-codelss  (GPS)

       I = I channel     (GPS, Galileo, QZSS, BDS)
       L = L channel     (L2C GPS, QZSS)
           P channel     (GPS. QZSS)
       M = M code-based  (GPS)
       N = Codeless      (GPS)
       P = P code-based  (GPS, GLONASS)
           Pilot channel (BDS)

       Q = Q channel     (GPS, Galileo, QZSS, BDS)
       S = D channel     (GPS, Galileo, QZSS, BDS)
           M channel     (L2C GPS, QZSS)

       W = Based on Z-tracking (GPS)
       X = B+C channels  (Galileo, IRNSS)
           I+Q channels  (GPS, IRNSS)
           M+L channels  (GPS, QZSS)
           D+P channels  (GPS, QZSS, BDS)

       Y = Y code based  (GPS)
       Z = A+B+C channels(Galileo)
           D+P channels  (BDS)
    -------------------------------------------------------------------------------------------------------------------------
    """
    eof         = 0
    success     = 1
    antDelta    = []
    timeSystem  = ''
    tFirstObs   = []
    tLastObs    = np.nan
    tInterval   = np.nan
    rinexProgr  = np.nan
    rinexDate   = np.nan
    obsCodes    = {}
    GNSSsystems = {}
    gnssType    = ""
    markerName  = ""
    numHeaderLines  = 0
    clockOffsetsON  = 0
    numGNSSsystems  = 0
    leapSec         = np.nan
    numOfObsCodes   = []
    rinexHeader     = {}
    approxPosition  = [0, 0, 0]
    obsCodeIndex = {}
    rinexVersion = np.nan
    recType = np.nan
    GLO_Slot2ChannelMap = np.nan

    ## -------Testing input arguments
    # Test if filename is valid format
    if type(filename) != str:
        print('INPUT ERROR(rinexReadsObsHeader211): The input argument filename is of type %s.\n Must be of type string or char' %(type(filename)))
        success = 0
        fid     = 0
        return success


    ## -- Open rinex observation file
    fid = open(filename,'r')
    if os.stat(filename).st_size == 0:
        raise ValueError('ERROR: This file seems to be empty')

    while 1: # Gobbling the header
        numHeaderLines = numHeaderLines + 1
        line = fid.readline().rstrip()

        if 'END OF HEADER' in line:
            break

        if numHeaderLines == 1: # if first line of header
            # store rinex version
            rinexVersion = line[0:9].strip()
            # store rinex type, ex. "N" or "O"
            rinexType = line[20]
            # if rinex file is not an observation file
            if rinexType != 'O':  # Rinex file is oservation file
                print('ERROR(rinexReadObsFileHeader211): the file is not a RINEX observations data file!')
                success = 0
                fid.close()
                return

            ## -- Check gnss type  ## Changend indent here 09.12.2022 (was apart of the if test above earlier, and thats wrong)
            gnssType = line[40] # reads the GNSS system type
            if gnssType not in [' ', 'G', 'R', 'C', 'E', 'M' ]:
                if gnssType in ['J', 'I', 'S']:
                    print('ERROR(rinexReadObsFileHeader211): This software is meant for reading GNSS data only.\
                           %s is an invalid satellite system type.' %(gnssType))
                else:
                    print('ERROR(rinexReadObsFileHeader211): %s is an unrecognized satellite system type.' %(gnssType))

                success = 0
                fid.close()
            ## -- If no system type, set G
            if gnssType == ' ':
                gnssType = 'G'


        if 'PGM / RUN BY / DATE' in line:
            rinexProgr = line[0:20] # rinex program
            rinexDate = line[40:60] # rinex date

        if 'MARKER NAME' in line:
            markerName = line.strip() # markername

        ## if no marker name, "MARKER" is read, so set to blank
        if 'Marker' in markerName:
            markerName = ''

        if 'ANTENNA: DELTA H/E/N' in line:
            for k in np.arange(0,3):
                line_ = [el for el in line.split(" ") if el != ""]
                antDelta = [line_[0],line_[1],line_[2]]


        if '# / TYPES OF OBSERV' in line:
            pattern_ep1 = r"\s*(\d{2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2})\s+(\d{1,2}\.\d+)\s*"
            pattern_ep2 = r'\s+([A-Z]\d{2})' #pattern for next line in block header
            first_block_lines = find_match_in_file_FWD(filename, pattern_ep1)
            sat_lines = [line for line in first_block_lines if re.match(pattern_ep1, line) or re.match(pattern_ep2, line)]
            sat_lines = ''.join([line.strip() for line in sat_lines if re.match(pattern_ep1, line) or re.match(pattern_ep2, line)])[29:] #[29:] is remove the first part of the string
            GNSSsystems = list(set(re.findall(r"[a-zA-Z]", sat_lines)))
            GNSSsystems = {i+1: GNSSsystems[i] for i in range(len(GNSSsystems))} # make dict there index in list becomes key in dict
            numGNSSsystems = len(GNSSsystems)

            line_dum = []
            while '# / TYPES OF OBSERV' in line:
                line = line[0:60]     # deletes 'SYS / # / OBS TYPES'
                line_ = [el for el in line.split(" ") if el != ""]
                line_dum.extend(line_)
                line = fid.readline().rstrip()
            line_ = line_dum
            nObs = int(line_.pop(0)) # assingning nObs to variable and removing it from the list
            undesiredobsCodeIndex = []
            desiredObsCodeIndex = []
            ## is Sys amoung desired GNSS systems
            GNSSSystemObsCodes = {}  # Reset cell of obsCodes for this GNSS system
            obsCode_list = []
            for k in np.arange(0,nObs):
                obsCode = line_.pop(0)
                # Checking if obsCode is valid
                if len(obsCode) != 2 or obsCode[0] not in ['C', 'L', 'D','S', 'I', 'X','P'] or  \
                          obsCode[1] not in ['0', '1', '2','3', '4', '5', '6', '7', '8', '9']:
                    print('ERROR (rinexReadsObsHeader211):  obsCode %s is a not a standard RINEX 2.11 observation type!' %(obsCode))

                ## is obsCode amoung desired obscodes and frequency bands
                if includeAllObsCodes or obsCode[0] in desiredObsCodes and int(obsCode[1]) in desiredObsBands:
                     ## store obsCode if amoung desire obsCodes
                    obsCode_list.append(obsCode)
                    for sys in GNSSsystems.values():
                        GNSSSystemObsCodes[sys] =  obsCode_list
                        GNSSSystemObsCodes[sys] =  obsCode_list

                    desiredObsCodeIndex.append(k)
                else:
                    # store index of discareded obsCode
                    undesiredobsCodeIndex.append(k)

            for sys_indx in GNSSsystems.keys():
                sys = GNSSsystems[sys_indx]
                numOfObsCodes.append(len(GNSSSystemObsCodes[sys]))
                obsCodes[sys_indx] = GNSSSystemObsCodes
                obsCodes[sys_indx] = GNSSSystemObsCodes

            obsCodeIndex[numGNSSsystems] = desiredObsCodeIndex # Store indices of desired obsCodes

        if 'PRN / # OF OBS' in line:
            system_list = []
            while 'PRN / # OF OBS' in line:
                line = line[0:60]     # deletes 'SYS / # / OBS TYPES'
                line_ = [el for el in line.split(" ") if el != ""]
                Sys = line_.pop(0)[0] # assingning nObs to variable and removing it from the list
                if Sys not in ["G","R","E","C"]: # added this line 29.01.2023 to fix bug where Only one system and several lines with Obscodes in rinex file
                    continue
                if (includeAllGNSSsystems and Sys in ["G", "R", "E", "C"] or Sys in desiredGNSSsystems):
                    # numGNSSsystems  = numGNSSsystems + 1 # increment number of GNSS systems
                    system_list  = [Sys] + [s for s in system_list if s not in [Sys]]
                    numGNSSsystems  = len(system_list)

                    if Sys not in GNSSsystems.values():
                        GNSSsystems[numGNSSsystems] = str(Sys) # Store current GNSS system
                    else:
                        pass
                    # GNSSsystems[numGNSSsystems] = str(Sys) # Store current GNSS system
                    GNSSSystemObsCodes[Sys] =  obsCode_list
                    numOfObsCodes.append(len(GNSSSystemObsCodes[Sys]))
                    obsCodes[numGNSSsystems] = GNSSSystemObsCodes
                    obsCodeIndex[numGNSSsystems] = desiredObsCodeIndex # Store indices of desired obsCodes
                    GNSSSystemObsCodes[Sys] =  obsCode_list
                numHeaderLines = numHeaderLines + 1
                line = fid.readline().rstrip()
                line_ = [el for el in line.split(" ") if el != ""]


        if 'TIME OF FIRST OBS' in line:
            line = line[0:60]     #  deletes 'TIME OF FIRST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            for k in np.arange(0,6):
                tok = line_.pop(0)  # finds the substrings containing the components of the time of the first observation
                                      #(YYYY; MM; DD; hh; mm; ss.sssssss) and specifies
                                      # the Time System used in the
                                      # observations file (GPST, GLOT or
                                      # GALT)
                if k ==0:
                    yyyy = int(tok)
                elif k ==1:
                    mm = int(tok)
                elif k ==2:
                    dd = int(tok)
                elif k ==3:
                    hh = int(tok)
                elif k ==4:
                    mnt = int(tok)
                elif k ==5:
                    ss = float(tok)


            tFirstObs = np.array([[yyyy],[mm],[dd],[hh],[mnt],[ss]])

            # Get Time system
            aux = line_.pop(0)
            # aux = strtok(line);
            if aux == 'GPS':
                timeSystem = 'GPS'
            elif aux == 'GLO':
                timeSystem = 'GLO'
            elif aux == 'GAL':
                timeSystem = 'GAL'
            elif aux == 'BDT':
                timeSystem = 'BDT'

            else:
                if gnssType == 'G':
                    timeSystem = 'GPST'
                elif gnssType == 'R':
                    timeSystem = 'GLOT'
                elif gnssType == 'E':
                    timeSystem = 'GALT'
                elif gnssType == 'C':
                    timeSystem = 'BDT'
                else:
                    print('CRITICAL ERROR (rinexReadsObsHeader211):\n' \
                                       'The Time System of the RINEX observations file '\
                                       'isn''t correctly specified!\n')
                    success = 0
                    fid.close()


        if 'TIME OF LAST OBS' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            for k in np.arange(0,6):
                tok = line_.pop(0)
                if k ==0:
                    yyyy = int(tok)
                elif k ==1:
                    mm = int(tok)
                elif k ==2:
                    dd = int(tok)
                elif k ==3:
                    hh = int(tok)
                elif k ==4:
                    mnt = int(tok)
                elif k ==5:
                    ss = float(tok)

            tLastObs = np.array([[yyyy],[mm],[dd],[hh],[mnt],[ss]])

        if 'INTERVAL' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            tInterval = float(line_.pop(0))


          ## -- This is an optional record!
          # if 'RCV CLOCK OFFS APPL' in line:
          #     if (strtok(line)=='0'):
          #         clockOffsetsON = 0;
          #     elif (strtok(line)=='1'):
          #         clockOffsetsON = 1;
          #     else:
          #         success = 0;
          #         print('ERROR (rinexReadsObsHeader211): unrecognized receiver clock offsets flag!')
          #         fid.close()


           ## This is an optional record
        if 'LEAP SECONDS' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            leapSec = int(line_.pop(0))


           ## -- store approximate receiver position
        if 'APPROX POSITION XYZ' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            approxPosition = np.array([[float(line_[0])],[float(line_[1])],[float(line_[2])]])


         ## GLOANSS SLOTS
        if 'GLONASS SLOT / FRQ #' in line:
            line = line[0:60]     #  deletes 'TIME OF LAST OBS'
            line_ = [el for el in line.split(" ") if el != ""]
            nGLOSat = int(line_.pop(0))
            slotNumbers = np.array([])
            channels = np.array([])
            for k in np.arange(0,nGLOSat):
                slotNumber = line_.pop(0)[1::]
                channel = int(line_.pop(0))
                slotNumbers = np.append(slotNumbers,slotNumber)
                channels = np.append(channels,channel)

                if np.mod(k+1, 8) == 0 and k+1 != 24:
                    # line = fgetl(fid); # end of line is reached so read next line
                    line = fid.readline().rstrip()
                    numHeaderLines = numHeaderLines + 1
                    line = line[0:60]     #  deletes 'TIME OF LAST OBS'
                    line_ = [el for el in line.split(" ") if el != ""]
                elif np.mod(k+1, 8) == 0 and k+1 == 24:
                    break

            GLO_Slot2ChannelMap = dict(zip(slotNumbers.astype(int),channels.astype(int)))

        if 'REC # / TYPE / VERS' in line:
            recType = line[20:40]


     # End of Gobbling Header Loop
    for k in np.arange(1,numGNSSsystems+1):
        # Give info if any of GNSS systems had zero of desired obscodes.
        if numOfObsCodes[k-1] == 0 or sum(tFirstObs) == 0:
            if GNSSsystems[k] == 'G':
                print('INFO: (rinexReadsObsHeader211)\nNone of the GPS satellites had any of the desired obsCodes\n\n')
            elif GNSSsystems[k] == 'R':
                print('INFO: (rinexReadsObsHeader211)\nNone of the GLONASS satellites had any of the desired obsCodes\n\n')

    ## store rinex header info
    rinexHeader['rinexVersion'] =rinexVersion
    rinexHeader['rinexType'] = rinexType
    rinexHeader['gnssType'] =gnssType
    rinexHeader['rinexProgr'] =rinexProgr
    rinexHeader['rinexDate'] =rinexDate


    print('INFO(rinexReadObsFileHeader211): Rinex header has been read')

    return success, rinexVersion, gnssType, markerName, recType, antDelta, GNSSsystems, numOfObsCodes, \
    obsCodes, obsCodeIndex,tFirstObs, tLastObs, tInterval,timeSystem, numHeaderLines, clockOffsetsON, \
    rinexProgr, rinexDate,leapSec, approxPosition, GLO_Slot2ChannelMap, eof, fid


def rinexReadObsBlock211(fid, numSV, nObsCodes, GNSSsystems, obsCodeIndex, readSS, readLLI, SVlist):
    """
    Reads all the observations from a RINEX observation block.

    Positioned at the beginning of the line immediately after the header of the
    observations block, reads all the observations in this block of a RINEX
    observations file. This function is meant to be used after using function
    rinexReadObsFileHeader211

    Based in the work of Antonio Pestana, rinexReadObsBlock211, March 2015
    --------------------------------------------------------------------------------------------------------------------------
    INPUTS:
    -------

    fid:                  Matlab file identifier of a Rinex observations text file

    numSV:                total number of satellites with observations in
                          current observation block, integer

    numOfObsCodes:        column vector containing number of observation
                          types for each GNSS system. Order is the same as
                          GNSSsystems

    GNSSsystems:          cell array containing codes of GNSS systems included
                          in RINEX observationfile. Elements are strings.
                          ex. "G" or "E"

    obsCodeIndex:         cell with one cell element for each GNSS system.
                          Order is the same as GNSSsystems. Each cell element
                          contains an array of indices. These indices
                          indicate the observation types that should be
                          read for each GNSS system. ex. If one index for
                          GPS is 1 then the first observation type for GPS
                          should be read.

    readSS:                   Boolean, 0 or 1.
                              1 = read "Signal Strength" Indicators
                              0 = do not read "Signal Strength" Indicators

    readLLI:                  Boolean, 0 or 1.
                              1 = read "Loss-Of-Lock Indicators"
                              0 = do not read "Loss-Of-Lock Indicators"
    --------------------------------------------------------------------------------------------------------------------------
    OUTPUTS:
    --------

    success:               Boolean. 1 if the function seems to be successful,
                          0 otherwise

    Obs:                  matrix [numSV x max_nObs] that stores all
                          observations of this observation block. max_nObs
                          is the highest number of observation codes that
                          any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    SVlist:               column cell [numSV x 1] that conatins the
                          identification code of each line of observation
                          block. ex. "G21". numSV is total number of
                          satellites minus amount of satellites removed.

    numSV:                numSV, unlike the input of same name, is the total
                          number of satellites minus amount of satellites
                          removed.

    LLI:                  matrix [numSV x max_nObs] that stores all
                          "loss-of-lock" indicators of this observation block.
                          max_nObs is the highest number of observation codes
                          that any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    SS:                   matrix [numSV x max_nObs] that stores all
                          "signal strength" indicators of this observation block.
                          max_nObs is the highest number of observation codes
                          that any of the GNSS systems have. Which observation
                          types that are associated with what collumn will
                          vary between GNSS systems. SVlist will give
                          overview of what GNSS system each row is connected
                          to

    eof:                  end-of-file flag; 1 if end-of-file was reached,
                          0 otherwise
    --------------------------------------------------------------------------------------------------------------------------
    """
    ## Initialize variables in case of input error
    success                     = np.nan
    eof                         = np.nan
    max_n_obs_Types             = np.nan
    Obs                         = np.nan
    LLI                         = np.nan
    SS                          = np.nan
    removed_sat                 = np.nan
    desiredGNSSsystems          = np.nan
    obsCodeIndex = obsCodeIndex[list(obsCodeIndex.keys())[0]]

    ## Test type of numSV
    if type(numSV) != int:
        print(f'INPUT ERROR(rinexReadObsBlock211): The input argument numSV is of type {type(numSV)}.\n Must be of type double')
        success = 0
        return success

    nObsCodes = int(nObsCodes[0])
    success = 1
    eof     = 0

    # Highest number of obs codes of any GNSS system
    max_n_obs_Types = nObsCodes
    # Initialize variables
    Obs = np.empty([numSV, max_n_obs_Types])
    if nObsCodes > 5:
        # There are three lines
        factor = (nObsCodes // 5) + 1
        nLines = factor * numSV
    else:
        factor = 1
        nLines = factor*numSV

    if readLLI:
        LLI = np.empty([numSV, max_n_obs_Types])

    if readSS:
        SS  = np.empty([numSV, max_n_obs_Types])

    # number of satellites excluded so far
    removed_sat = 0
    desiredGNSSsystems = list(GNSSsystems.values())
    #There are cases where matching negative values occurs
    pattern = r"^\s*-?[0-9]+\.[0-9]+\s*"
    pattern2 = r'\s*-?\d+\.\d+'
    pattern4 = r'^\s*[A-Za-z]+[A-Za-z0-9]*(?:[\s]+[A-Za-z]+[A-Za-z0-9]*)*\s*$'
    ## -- Gobble up observation block
    for sat in np.arange(0,numSV):
        line = fid.readline().rstrip()
        while (not re.match(pattern,line) or not re.match(pattern2,line) or re.match(pattern4,line)) and line !="":
            sat_overview = line.split(" ")[-1]
            pattern3 = re.compile(r'[A-Z][0-9]{2}')
            sat_list = re.findall(pattern3, sat_overview)
            for s in sat_list:
                sys = s[0]
                if sys not in desiredGNSSsystems:
                    continue
                SVlist.append(s)
            line = fid.readline().rstrip()

        SV = SVlist[sat]
        if SV[0] not in desiredGNSSsystems:
            removed_sat +=1
        else:
            ## Index of current GNSS system
            GNSSsystemIndex = [i for i in GNSSsystems if GNSSsystems[i]==SV[0]][0]
            n_obs_current_system = nObsCodes
            nNew_line = 1 # counter to keep track in new line of same satellite
            for obs_num in np.arange(0, n_obs_current_system):
                # obsIndex = obsCodeIndex[GNSSsystemIndex][obs_num]
                obsIndex = obsCodeIndex[obs_num]
                # charPos = 4+(obsIndex)*16
                if obsIndex < 5:
                    charPos = 1+(obsIndex)*16
                else:
                    # There are three lines
                    charPos = 1 + (obsIndex % 5) * 16

                ## check that the current observation of the current GNSS system
                ## is not on the list of obs types to be excluded
                ## stringlength of next obs.
                obsLen = min(14, len(line) - charPos)
                # read next obs
                newObs = line[charPos:charPos+obsLen].strip()
                # If observation missing, set to 0
                if newObs != '':
                    newObs = float(newObs)
                else:
                    newObs = 0
               # Store new obs
                Obs[sat - removed_sat, obs_num] = newObs

                if readLLI:
                # read LLI of current obs (if present)
                    if charPos+13<len(line): ## endret til < (kun mindre enn) 13.11
                        newLLI = line[charPos+13] # loss of lock indicator ### endret fra 14 til 13 den 13.11.2022
                    else:
                        newLLI = ' '
                    if newLLI.isspace():
                        newLLI = -999
                    else:
                        newLLI = int(newLLI)
                     # Store LLI
                    LLI[sat - removed_sat, obs_num] = newLLI

                if readSS:
                    # read SS of current obs (if present)
                    if charPos+14<len(line): ## endret til < (kun mindre enn) 13.11
                        newSS = line[charPos+14] # signal strength endret fra 15 til 14 den 13.11.2022
                    else:
                        newSS = ' '
                    # if no SS set to -999
                    if newSS.isspace():
                        newSS = -999
                    else:
                        newSS = int(newSS)

                    ## Store SS
                    SS[sat - removed_sat, obs_num]  = newSS

                    # if np.mod(obs_num+1, 5) == 0 and nObsCodes > 5 and factor*sat < factor*numSV:
                    # Matches the case where three lines are read
                    if np.mod(obs_num+1, 5) == 0 and nObsCodes>5 and nNew_line <= factor:
                        nNew_line += 1
                        line = fid.readline().rstrip()


    ## -- Update number og satellites after satellites have been excluded
    numSV = numSV - removed_sat
    ## --Remove empty arrays
    idx_keep = len(Obs) -1 -removed_sat + 1 # removing sats
    Obs = Obs[:idx_keep,:]
    return success, Obs, SVlist, numSV, LLI, SS, eof


def rinexReadObsBlockHead211(fid):
    """
    Reads the metadata in the head of a RINEX 3.xx observations block, NOT
    the header of the file.

    ATTENTION: Ignores all data in blocks with event flags with numbers
    greater than 1!!!

    Positioned in a RINEX V.2 GNSS observations text file at the beginning
    of an observation block. In rinex 3.xx the line starts with '> '

    --------------------------------------------------------------------------------------------------------------------------
     INPUTS

     fid:              Matlab identifier of an open RINEX V.2 GNSS
                       observations text file positioned at the beginning
                       of an observation block.
    --------------------------------------------------------------------------------------------------------------------------
     OUTPUTS

     success:          1 if function performs successfully, 0 otherwise

     epochflag:        Rinex observations epoch flag, as follows:
                           0: OK
                           1: power failure between previous and current epoch
                       From now on the "event flags":
                           2: start moving antenna
                           3: new site occupation
                           4: header information follows
                           5: external event (epoch is significant)

     clockOffset:          value of the receiver clock offset. If not present
                           in the metadata of the observations block
                           (it's optional RINEX V.2 data)it is assumed to be
                           zero. If not zero implies that epoch, code, and
                           phase data have been corrected by applying
                           realtime-derived receiver clock offset

     date:                 time stamp of the observations block. Six-elements column-vector
                           as follows:
                               year: four-digits year (eg: 1959)
                               month: integers 1..12
                               day: integers 1..31
                               hour: integers 0..24
                               minute: integers 0..60
                               second: reals 0..60

     SVlist:               Dictionary with the system and PRN number for the current epoch in obsfile. Ex {G:[3,4,20,24,28]
                                                                                                           R:[1,3,4,9,12]}

     numSV:                number of satellites with observations in with
                           observations. This will include all satellite
                           systems.
    --------------------------------------------------------------------------------------------------------------------------

    """
    ## -- Initialize variables
    success = 1
    eof     = 0
    date    = [0,0,0,0,0,0]
    numSV   = 0
    epochflag = 0
    clockOffset = 0
    noFlag = 1
    SVlist = np.nan
    line = fid.readline().rstrip()
    if not line:
        eof = 1
        print('\nINFO(rinexReadObsBlockHead211): End of observations text file reached')
        return success, epochflag, clockOffset, date, numSV,SVlist, eof

    epochflag   = int(line[28])
    # Flags 2-5 are special events (moved antenna, header follows, etc.)
    # The number of special records to skip is at columns 30-32.
    while epochflag > 1:
        noFlag = 0
        num_special = int(line[29:32])
        for _ in range(num_special):
            fid.readline()
        msg = 'WARNING(rinexReadsObsBlockHead211): Observations event flag %d encountered. %d special record(s) skipped.' % (epochflag, num_special)
        print(msg)
        # Read the next epoch header line
        line = fid.readline().rstrip()
        if not line:
            eof = 1
            print('\nINFO(rinexReadObsBlockHead211): End of observations text file reached')
            return success, epochflag, clockOffset, date, numSV, SVlist, eof
        epochflag = int(line[28])

    # Gets the number of used satellites in obs epoch
    numSV = int(line[30:32])
    # Gets the receiver clock offset. This is optional data!
    # clockOffset = 0
    # if len(line) == 56:
    #     clockOffset = float(line[41:56])

    ## -- Reads the time stamp of the observations block (6 numerical values)
    date = [el for el in line[1::].split(" ") if el != ""]
    date = date[:6]
    date = [float(el) for el in date]

    if noFlag == 0:
        # RINEX 2 uses 2-digit years; convert for display
        yr = int(date[0])
        yr4 = yr + 2000 if yr < 80 else yr + 1900
        print('Epoch date = %.4d %.2d %.2d %.2d:%.2d:%6.4f' % (yr4,date[1],date[2],date[3],date[4],date[5]))


    SVlist = re.findall(r'[A-Z][0-9]{2}', line)

    return success, epochflag, clockOffset, date, numSV, SVlist, eof
