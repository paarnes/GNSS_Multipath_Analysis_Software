try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    from importlib_metadata import version, PackageNotFoundError  # type: ignore

try:
    __version__ = version("gnssmultipath")
except PackageNotFoundError:
    __version__ = "unknown"

from .GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
from .readRinexObs import readRinexObs
from .RinexNav import Rinex_v3_Reader
from .RinexNav import Rinex_v2_Reader
from .PickleHandler import PickleHandler
from .read_SP3Nav import readSP3Nav
from .Geodetic_functions import *
from .SatelliteEphemerisToECEF import SatelliteEphemerisToECEF
from .GNSSPositionEstimator import GNSSPositionEstimator
from .SP3PositionEstimator import SP3PositionEstimator