from typing import Literal, Dict, Tuple, Optional
import numpy as np
from gnssmultipath.BroadNavPositionEstimator import BroadNavPositionEstimator
from gnssmultipath.SP3PositionEstimator import SP3PositionEstimator
from pyproj import Transformer


class GNSSPositionEstimator:
    """
    A class to estimate GNSS receiver positions using either navigation data or precise SP3 files.


    Parameters:
    ----------
    rinex_obs_file : str
        Path to the RINEX observation file.
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
                - If an invalid CRS is provided, the transformation defaults to 'ECEF'.
                - Use https://epsg.org/home.html to find the EPSG code for the desired CRS.

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

    """

    def __init__(
        self,
        rinex_obs_file: str,
        desired_time: np.ndarray,
        rinex_nav_file: str = None,
        sp3_file: str = None,
        desired_system: Literal["G", "E", "R"] = "G",
        x_rec_approx: float = 0.0,
        y_rec_approx: float = 0.0,
        z_rec_approx: float = 0.0,
        elevation_cut_off_angle: Optional[int] = 10,
        crs: str = "EPSG:4978",
    ):
        if rinex_nav_file:
            self.GNSSPos = BroadNavPositionEstimator(
                rinex_obs_file=rinex_obs_file,
                rinex_nav_file=rinex_nav_file,
                desired_system=desired_system,
                desired_time=desired_time,
                x_rec_approx=x_rec_approx,
                y_rec_approx=y_rec_approx,
                z_rec_approx=z_rec_approx,
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

        try:
            transformer = Transformer.from_crs("EPSG:4978", self.crs, always_xy=True)
            transformed_coords = transformer.transform(coords[0], coords[1], coords[2])
            return np.array([transformed_coords[0], transformed_coords[1], transformed_coords[2], coords[3]])
        except Exception as e:
            print(f"Error: Transformation not possible! Returning receiver coordinates in ECEF (EPSG:4978). Error message: {e}")
            return coords


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





