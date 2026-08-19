try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    from importlib_metadata import version, PackageNotFoundError  # type: ignore

try:
    __version__ = version("gnssmultipath")
except PackageNotFoundError:
    # Fallback: read version from pyproject.toml when not installed as a package
    import pathlib as _pathlib
    _pyproject = _pathlib.Path(__file__).resolve().parents[2] / "pyproject.toml"
    if _pyproject.is_file():
        import re as _re
        _match = _re.search(r'^version\s*=\s*"([^"]+)"', _pyproject.read_text(), _re.MULTILINE)
        __version__ = _match.group(1) if _match else "unknown"
    else:
        __version__ = "unknown"

from .GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
from .signalAnalysis import SignalAnalyzer, SignalStats, SlipPeriods
from .readers.readRinexObs import readRinexObs, RinexObsData
from .readers.GNSSObservationData import (
    GNSSObservationData, SystemObservations, SatelliteObservations,
    EpochObservations,
)
from .readers.ObsCode import ObsCode
from .readers.RinexNav import RinexNav, RinexNavData
from .utils.PickleHandler import PickleHandler
from .readers.read_SP3Nav import readSP3Nav
from .readers.SP3Reader import SP3Reader
from .constants import SPEED_OF_LIGHT, carrier_frequency, wavelength
from .Geodetic_functions import *
from .SatelliteEphemerisToECEF import SatelliteEphemerisToECEF
from .GNSSPositionEstimator import GNSSPositionEstimator
from .SP3PositionEstimator import SP3PositionEstimator
from .SP3Interpolator import SP3Interpolator
from .PreciseSatCoords import PreciseSatCoords
from .utils import CDDISDownloader