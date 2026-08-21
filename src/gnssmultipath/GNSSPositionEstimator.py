from typing import Literal, Dict, Tuple, Optional, Union
import numpy as np
import pandas as pd
from gnssmultipath.BroadNavPositionEstimator import BroadNavPositionEstimator
from gnssmultipath.SP3PositionEstimator import SP3PositionEstimator
from gnssmultipath.SatelliteEphemerisToECEF import SatelliteEphemerisToECEF
from gnssmultipath.readers.readRinexObs import RinexObsData
from pyproj import Transformer


class GNSSPositionEstimator:
    """
    A class to estimate GNSS receiver positions using either navigation data or precise SP3 files.


    Parameters:
    ----------
    rinex_obs_file : str
        Path to the RINEX observation file. Can be omitted when ``rinex_data``
        is supplied.
    desired_time : np.ndarray
        Desired observation time as an array [year, month, day, hour, minute, second].
    rinex_nav_file : str, optional
        Path to the RINEX navigation file (default is None).
    sp3_file : str, optional
        Path to the SP3 precise ephemeris file (default is None).
    desired_system : Literal['G', 'E', 'R'], optional
        Desired GNSS system: 'G' (GPS), 'E' (Galileo), 'R' (GLONASS) (default is 'G').
    x_rec_approx : float, optional
        Approximate X-coordinate of the receiver (default is 0).
    y_rec_approx : float, optional
        Approximate Y-coordinate of the receiver (default is 0).
    z_rec_approx : float, optional
        Approximate Z-coordinate of the receiver (default is 0).
    elevation_cut_off_angle : int, optional
        Elevation cut-off angle for satellites (default is 10).
    crs : str, optional
        Coordinate Reference System (CRS) to define the coordinate format of the receiver's position.
                - Default is EPSG:4978 which is WGS84 Cartesian 3D coordinates(ECEF).
                - To get longitude, latitude, and altitude (LLA) coordinates, use EPSG:4326
                - The default needs no transformation, so PROJ is only required when
                  another CRS is requested.
                - If the transformation fails, a RuntimeError is raised rather than
                  silently returning ECEF coordinates.
                - Use https://epsg.org/home.html to find the EPSG code for the desired CRS.
    navdata : SatelliteEphemerisToECEF, optional
        Already converted broadcast ephemerides. Supplying this avoids reading the
        navigation file again.
    rinex_data : RinexObsData, optional
        Already parsed RINEX observation data. When provided, its fields are
        used instead of reading ``rinex_obs_file``.
    sp3_metadata_dict : dict, optional
        SP3 metadata, required when ``sp3_file`` is an already read pandas DataFrame.

    Example:
    -------
    .. code-block:: python
            rinObs = "your_rinex_obs_file"
            sp3 = "your_sp3_file"
            desired_time = np.array([2022, 1, 1, 0, 0, 30.0000000])
            desired_system = "G"
            gnsspos, stats = GNSSPositionEstimator(rinObs,
                                            sp3_file = sp3,
                                            desired_time = desired_time,
                                            desired_system = desired_system,
                                            crs = "EPSG:32632", # WGS84 UTM Zone 32N
                                            ).estimate_postiion()


    Reusing data that is already read in:

    .. code-block:: python
            rinex = readRinexObs(rinObs)
            navdata = SatelliteEphemerisToECEF(rinNav, 0, 0, 0)
            gnsspos, stats = GNSSPositionEstimator(desired_time=desired_time,
                                            navdata=navdata,
                                            rinex_data=rinex,
                                            ).estimate_position()

    """

    def __init__(
        self,
        rinex_obs_file: Optional[str] = None,
        desired_time: np.ndarray = None,
        rinex_nav_file: str = None,
        sp3_file: Union[str, "pd.DataFrame"] = None,
        desired_system: Literal["G", "E", "R"] = "G",
        x_rec_approx: float = 0.0,
        y_rec_approx: float = 0.0,
        z_rec_approx: float = 0.0,
        elevation_cut_off_angle: Optional[int] = 10,
        crs: str = "EPSG:4978",
        navdata: Optional[SatelliteEphemerisToECEF] = None,
        sp3_metadata_dict: Optional[dict] = None,
        rinex_data: Optional[RinexObsData] = None,
    ):
        if rinex_obs_file is not None and rinex_data is not None:
            raise ValueError(
                "Provide either rinex_obs_file or rinex_data, not both"
            )

        GNSS_obs = time_epochs = GNSSsystems = obsCodes = None
        if rinex_data is not None:
            if not isinstance(rinex_data, RinexObsData):
                raise TypeError("rinex_data must be a RinexObsData instance")
            rinex_obs_file = None
            GNSS_obs = rinex_data.GNSS_obs
            time_epochs = rinex_data.time_epochs
            GNSSsystems = rinex_data.GNSSsystems
            obsCodes = rinex_data.obsCodes

        if rinex_nav_file or navdata is not None:
            self.GNSSPos = BroadNavPositionEstimator(
                rinex_obs_file=rinex_obs_file,
                rinex_nav_file=rinex_nav_file,
                navdata=navdata,
                desired_system=desired_system,
                desired_time=desired_time,
                x_rec_approx=x_rec_approx,
                y_rec_approx=y_rec_approx,
                z_rec_approx=z_rec_approx,
                GNSS_obs=GNSS_obs,
                time_epochs=time_epochs,
                GNSSsystems=GNSSsystems,
                obsCodes=obsCodes,
                elevation_cut_off_angle = elevation_cut_off_angle,
            )
        else:
            self.GNSSPos = SP3PositionEstimator(
                sp3_data=sp3_file,
                rinex_obs_file=rinex_obs_file,
                desired_time=desired_time,
                desired_system=desired_system,
                x_rec_approx=x_rec_approx,
                y_rec_approx=y_rec_approx,
                z_rec_approx=z_rec_approx,
                GNSS_obs=GNSS_obs,
                time_epochs=time_epochs,
                GNSSsystems=GNSSsystems,
                obsCodes=obsCodes,
                sp3_metadata_dict=sp3_metadata_dict,
                elevation_cut_off_angle = elevation_cut_off_angle,
            )

        self.crs = crs

    def __transform_coordinates(self, coords: np.ndarray) -> np.ndarray:
        """
        Transform the estimated coordinates to the desired coordinate reference system.

        Parameters
        ----------
        coords : np.ndarray
            Estimated coordinates as a numpy array [X, Y, Z, dT].

        Returns
        -------
        transformed_coords : np.ndarray
            Transformed coordinates to the desired coordinate reference system.
        """

        # The position is already estimated in EPSG:4978, so the default case
        # needs no transformation - and therefore no working PROJ installation.
        if str(self.crs).strip().upper().replace(" ", "") in ("EPSG:4978", "4978"):
            return coords

        try:
            transformer = Transformer.from_crs("EPSG:4978", self.crs, always_xy=True)
            transformed_coords = transformer.transform(coords[0], coords[1], coords[2])
        except Exception as e:
            # Returning ECEF here would hand back numbers in a completely
            # different frame than the caller asked for, so fail instead.
            raise RuntimeError(
                f"Could not transform the estimated position from EPSG:4978 to '{self.crs}'. "
                f"Check that the CRS code is valid and that PROJ is installed correctly. "
                f"Original error: {e}"
            ) from e

        if not all(np.isfinite(c) for c in transformed_coords[:3]):
            raise RuntimeError(
                f"Transformation from EPSG:4978 to '{self.crs}' produced a non-finite "
                'result. Check that the CRS is appropriate for the receiver location.'
            )

        return np.array([transformed_coords[0], transformed_coords[1], transformed_coords[2], coords[3]])


    def estimate_position(self) -> Tuple[np.ndarray, Dict]:
        """
        Estimate the GNSS receiver position and associated statistical data.

        Returns
        -------
        estimated_position : np.ndarray
            Estimated receiver position as a numpy array [X, Y, Z, dT].
        stats : Dict
            Dictionary containing statistical information, including:
                - Residuals
                - Sum of squared errors (SSE)
                - Standard deviation of unit weight (S0)
                - Cofactor and covariance matrices
                - Dilution of Precision (DOP) values
                - Standard deviations of parameters
        """
        coords, stats = self.GNSSPos.estimate_position()

        if self.crs != "ECEF":
            coords = self.__transform_coordinates(coords)
        return coords, stats





if __name__ == "__main__":
    pass





