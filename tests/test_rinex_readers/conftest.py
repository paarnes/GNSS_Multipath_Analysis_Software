"""
Shared fixtures for RINEX reader tests.
"""
import sys
import os
import pytest

# Setup project paths
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(project_path, 'src'))

TESTDATA_DIR = os.path.join(project_path, "TestData")
OBS_DIR = os.path.join(TESTDATA_DIR, "ObservationFiles")
NAV_DIR = os.path.join(TESTDATA_DIR, "NavigationFiles")
SP3_DIR = os.path.join(TESTDATA_DIR, "SP3")


@pytest.fixture(scope="module")
def rinex304_obs_file():
    return os.path.join(OBS_DIR, "v3", "OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx")


@pytest.fixture(scope="module")
def rinex304_obs_full_file():
    return os.path.join(OBS_DIR, "v3", "OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx")


@pytest.fixture(scope="module")
def rinex305_highrate_obs_file():
    return os.path.join(
        OBS_DIR, "v3", "BLUF00NZL_R_20251881200_02H_10Z_MO_first1min.rnx"
    )


@pytest.fixture(scope="module")
def rinex211_obs_file():
    return os.path.join(OBS_DIR, "v2", "gmgd31000_v2_11.20o")


@pytest.fixture(scope="module")
def rinex211_obs_file_2():
    return os.path.join(OBS_DIR, "v2", "p0803430_v211.24o")


@pytest.fixture(scope="module")
def rinex210_obs_file_event_flags():
    """RINEX 2.10 file with 10 obs codes per system (a multiple of 5) plus
    event-flag (2 / 5) special-record blocks. Truncated copy of an OPEC
    converter output that previously triggered an off-by-one continuation
    line read in ``rinexReadObsBlock211`` and crashed the v2 parser with
    ``ValueError: invalid literal for int() with base 10: '  '``.
    """
    return os.path.join(OBS_DIR, "v2", "OPEC0010_truncated.22o")


@pytest.fixture(scope="module")
def nav_gps_file():
    return os.path.join(NAV_DIR, "v3", "OPEC00NOR_S_20220010000_01D_GN.rnx")


@pytest.fixture(scope="module")
def nav_galileo_file():
    return os.path.join(NAV_DIR, "v3", "OPEC00NOR_S_20220010000_01D_EN.rnx")


@pytest.fixture(scope="module")
def nav_glonass_file():
    return os.path.join(NAV_DIR, "v3", "OPEC00NOR_S_20220010000_01D_RN.rnx")


@pytest.fixture(scope="module")
def nav_beidou_file():
    return os.path.join(NAV_DIR, "v3", "OPEC00NOR_S_20220010000_01D_CN.rnx")


@pytest.fixture(scope="module")
def nav_mixed_file():
    return os.path.join(NAV_DIR, "v3", "BRDC00IGS_R_20220010000_01D_MN.rnx")


@pytest.fixture(scope="module")
def nav_highrate_mixed_file():
    return os.path.join(NAV_DIR, "v3", "BRDC00IGS_R_20251880000_01D_MN.rnx")


@pytest.fixture(scope="module")
def nav_v4_mixed_file():
    return os.path.join(NAV_DIR, "v4", "BRD400DLR_S_20230710000_01D_MN_rin_v4.rnx")


@pytest.fixture(scope="module")
def nav_v2_gps_file():
    return os.path.join(NAV_DIR, "v2", "auto3430_v211.24n")


@pytest.fixture(scope="module")
def nav_v2_glonass_file():
    """RINEX v2.10 GLONASS broadcast nav (BRUX, 2024 DOY 001).

    Contains slots 1..25 -- exercises the v2 GLONASS reader path
    (4-line blocks, ``Rxx`` PRN labels, FCN extraction).
    """
    return os.path.join(NAV_DIR, "v2", "BRUX0010_v210.24g")


@pytest.fixture(scope="module")
def sp3_file_2022():
    return os.path.join(SP3_DIR, "Testfile_20220101.eph")


@pytest.fixture(scope="module")
def sp3_file_2020():
    return os.path.join(SP3_DIR, "Testfile_2020_11_5.SP3")


@pytest.fixture(scope="module")
def sp3_file_nmbus():
    return os.path.join(SP3_DIR, "NMBUS_2020 10 30.SP3")
