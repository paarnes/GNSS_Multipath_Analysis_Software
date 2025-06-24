"""
This module is using pytest to run tests on the software.
Each pull request to the master branch will trigger
a workflow that run these tests.

Made by: Per Helge Aarnes
E-mail: per.helge.aarnes@gmail.com

Example on how to run the tests:

To run all tests:
    pytest test_GNSS_MultipathAnalysis.py  -vv

To run a specific test:
    pytest test_GNSS_MultipathAnalysis.py::test_GNSS_MultipathAnalysis_sp3_file  -vv

To run a specific test in debug mode to see where it fails:
    pytest test_GNSS_MultipathAnalysis.py::test_GNSS_MultipathAnalysis_sp3_file  -vv --pdb  (type "quit" to exit the debug mode)

"""

import sys
import os
import pytest
from numpy.testing import assert_almost_equal

project_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(os.path.join(project_path,'src'))

from gnssmultipath import GNSS_MultipathAnalysis
from gnssmultipath.Geodetic_functions import ECEF2geodb, ECEF2enu, compute_azimut_elev
from gnssmultipath.PickleHandler import PickleHandler


test_data_ECEF2geodb = [
    (6378137.0, 6356752.314245, 2765120.7658, -4449250.0340, -3626405.4228, (-0.6086611088458492, -1.014732114220556, 41.7607967723161)),
    (6378137.0, 6356752.314245, -1270826.8700, 6242631.4460, 307792.4399, (0.048601283617660994, 1.771624423298406, 15.002176421694458)),
    (6378137.0, 6356752.314245, -4052052.7300, 4212835.9796, -2545104.5870, (-0.413121356617561, 2.3367431730495127, 603.2359722275287))
]

test_data_ECEF2enu = [
    (0.2838307690924083, -1.0738580938997349, 12893444.612051928, -9888519.684510132, 13138615.445688058, (6619718.534261429, 8457426.340380616, 17924793.375904225))
]


test_data_compute_azimut_elev = [(15813230.324051928, -15272264.751510132, 14913220.137688057, 2919785.712, -5383745.067, 1774604.692, (38.05066871295592, 59.07109574929687)) #(X,Y,Z,xm,ym,zm, (azimut,elevation (expected)))
]


@pytest.mark.parametrize("a, b, X, Y, Z, expected_output", test_data_ECEF2geodb)
def test_ECEF2geodb(a, b, X, Y, Z, expected_output):
    result = ECEF2geodb(a, b, X, Y, Z)
    assert result == pytest.approx(expected_output, abs=1e-6)


@pytest.mark.parametrize("lat, lon, dX, dY, dZ, expected_output", test_data_ECEF2enu)
def test_ECEF2enu(lat, lon, dX, dY, dZ, expected_output):
    result = ECEF2enu(lat,lon,dX,dY,dZ)
    assert result == pytest.approx(expected_output, abs=1e-3)


@pytest.mark.parametrize("X, Y, Z, xm, ym, zm, expected_output", test_data_compute_azimut_elev)
def test_compute_azimut_elev(X, Y, Z, xm, ym, zm, expected_output):
    result = compute_azimut_elev(X, Y, Z, xm, ym, zm)
    assert result == pytest.approx(expected_output, abs=1e-4)


os.chdir(os.path.join(project_path,'src'))

@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_GNSS_MultipathAnalysis_sp3_file():
    """
    Test the results from OPEC2022 and sp3 files
    """

    rinObs_file = "../TestData/ObservationFiles/NMBUS_SAMSUNG_S20.20o"
    sp3Nav_file = "../TestData/SP3/NMBUS_2020 10 30.SP3"
    expected_res = "../tests/analysisResults.pkl.zst"
    result = GNSS_MultipathAnalysis(rinObsFilename=rinObs_file, sp3NavFilename_1=sp3Nav_file,
                                    plotEstimates=False,
                                    write_results_to_csv=False,
                                    plot_polarplot=False)

    expected_result = PickleHandler.read_zstd_pickle(expected_res)
    # Compare the result with the expected result
    assert_almost_equal(expected_result['Sat_position']['G']['position']["2"],  result['Sat_position']['G']['position']["2"],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['G']['position']["4"],  result['Sat_position']['G']['position']["4"],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['G']['position']["30"],  result['Sat_position']['G']['position']["30"],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['G']['azimuth'],  result['Sat_position']['G']['azimuth'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['azimuth'],  result['Sat_position']['E']['azimuth'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['G']['elevation'],  result['Sat_position']['G']['elevation'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['C']['azimuth'],  result['Sat_position']['C']['azimuth'],decimal=4)
    assert_almost_equal(expected_result['GPS']['Band_1']['C1C']['rms_multipath_range1_averaged'], result['GPS']['Band_1']['C1C']['rms_multipath_range1_averaged'], decimal=4)
    assert_almost_equal(expected_result['GPS']['Band_1']['C1C']['elevation_weighted_average_rms_multipath_range1'], result['GPS']['Band_1']['C1C']['elevation_weighted_average_rms_multipath_range1'], decimal=4)

    assert expected_result['GPS']['Band_1']['C1C']['cycle_slip_distribution'] == result['GPS']['Band_1']['C1C']['cycle_slip_distribution']







@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_GNSS_MultipathAnalysis_broadcast_navfile():
    """
    Test the results from OPEC2022 and use of broadcast ephemerids
    """
    # os.chdir('src')
    rinObs_file =  "../TestData/ObservationFiles/OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.rnx"
    broadNav_file = "../TestData/NavigationFiles/BRDC00IGS_R_20220010000_01D_MN.rnx"
    expected_res = "../tests/analysisResults_OPEC00NOR_S_20220010000_01D_30S_MO_3.04_croped.pkl.zst"

    result = GNSS_MultipathAnalysis(rinObsFilename=rinObs_file, broadcastNav1=broadNav_file,
                                    plotEstimates=False,
                                    plot_polarplot=False,
                                    write_results_to_csv=False,
                                    nav_data_rate=120)

    expected_result = PickleHandler.read_zstd_pickle(expected_res)
    # Compare the result with the expected result
    assert_almost_equal(expected_result['Sat_position']['G']['position']['1'],  result['Sat_position']['G']['position']['1'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['G']['position']['8'],  result['Sat_position']['G']['position']['8'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['G']['position']['10'], result['Sat_position']['G']['position']['10'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['R']['position']['1'],  result['Sat_position']['R']['position']['1'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['R']['position']['7'],  result['Sat_position']['R']['position']['7'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['R']['position']['8'],  result['Sat_position']['R']['position']['8'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['position']['1'],  result['Sat_position']['E']['position']['1'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['position']['3'],  result['Sat_position']['E']['position']['3'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['position']['8'],  result['Sat_position']['E']['position']['8'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['C']['position']['5'],  result['Sat_position']['C']['position']['5'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['C']['position']['6'],  result['Sat_position']['C']['position']['6'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['C']['position']['9'],  result['Sat_position']['C']['position']['9'],decimal=4)

    assert_almost_equal(expected_result['Sat_position']['G']['azimuth'],  result['Sat_position']['G']['azimuth'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['R']['azimuth'],  result['Sat_position']['R']['azimuth'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['azimuth'],  result['Sat_position']['E']['azimuth'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['C']['azimuth'],  result['Sat_position']['C']['azimuth'],decimal=4)

    assert_almost_equal(expected_result['Sat_position']['G']['elevation'],  result['Sat_position']['G']['elevation'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['R']['elevation'],  result['Sat_position']['R']['elevation'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['elevation'],  result['Sat_position']['E']['elevation'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['C']['elevation'],  result['Sat_position']['C']['elevation'],decimal=4)

    assert_almost_equal(expected_result['GPS']['Band_1']['C1C']['rms_multipath_range1_averaged'], result['GPS']['Band_1']['C1C']['rms_multipath_range1_averaged'], decimal=4)
    assert_almost_equal(expected_result['GPS']['Band_1']['C1C']['elevation_weighted_average_rms_multipath_range1'], result['GPS']['Band_1']['C1C']['elevation_weighted_average_rms_multipath_range1'], decimal=4)
    assert_almost_equal(expected_result['GLONASS']['Band_1']['C1C']['rms_multipath_range1_averaged'], result['GLONASS']['Band_1']['C1C']['rms_multipath_range1_averaged'], decimal=4)
    assert_almost_equal(expected_result['GLONASS']['Band_1']['C1C']['elevation_weighted_average_rms_multipath_range1'], result['GLONASS']['Band_1']['C1C']['elevation_weighted_average_rms_multipath_range1'], decimal=4)
    assert_almost_equal(expected_result['Galileo']['Band_1']['C1X']['rms_multipath_range1_averaged'], result['Galileo']['Band_1']['C1X']['rms_multipath_range1_averaged'], decimal=4)
    assert_almost_equal(expected_result['Galileo']['Band_1']['C1X']['elevation_weighted_average_rms_multipath_range1'], result['Galileo']['Band_1']['C1X']['elevation_weighted_average_rms_multipath_range1'], decimal=4)
    assert_almost_equal(expected_result['BeiDou']['Band_2']['C2X']['rms_multipath_range1_averaged'], result['BeiDou']['Band_2']['C2X']['rms_multipath_range1_averaged'], decimal=4)
    assert_almost_equal(expected_result['BeiDou']['Band_2']['C2X']['elevation_weighted_average_rms_multipath_range1'], result['BeiDou']['Band_2']['C2X']['elevation_weighted_average_rms_multipath_range1'], decimal=4)

    assert_almost_equal(expected_result['GPS']['Band_1']['C1C']['nEstimates'], result['GPS']['Band_1']['C1C']['nEstimates'], decimal=4)
    assert_almost_equal(expected_result['GLONASS']['Band_1']['C1C']['nEstimates'], result['GLONASS']['Band_1']['C1C']['nEstimates'], decimal=4)
    assert_almost_equal(expected_result['Galileo']['Band_1']['C1X']['nEstimates'], result['Galileo']['Band_1']['C1X']['nEstimates'], decimal=4)
    assert_almost_equal(expected_result['BeiDou']['Band_2']['C2X']['nEstimates'], result['BeiDou']['Band_2']['C2X']['nEstimates'], decimal=4)

    assert expected_result['GPS']['Band_1']['C1C']['cycle_slip_distribution']     == result['GPS']['Band_1']['C1C']['cycle_slip_distribution']
    assert expected_result['GLONASS']['Band_1']['C1C']['cycle_slip_distribution'] == result['GLONASS']['Band_1']['C1C']['cycle_slip_distribution']
    assert expected_result['Galileo']['Band_1']['C1X']['cycle_slip_distribution'] == result['Galileo']['Band_1']['C1X']['cycle_slip_distribution']
    assert expected_result['BeiDou']['Band_2']['C2X']['cycle_slip_distribution']  == result['BeiDou']['Band_2']['C2X']['cycle_slip_distribution']




if __name__ == '__main__':
    pytest.main()
