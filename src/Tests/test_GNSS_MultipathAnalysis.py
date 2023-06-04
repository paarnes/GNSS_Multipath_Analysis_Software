import pytest
import pickle
import sys
import os
import numpy as np
from numpy.testing import assert_almost_equal

# os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\Multipath_analysis")
sys.path.append(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\Multipath_analysis")
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis
os.chdir(r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\Multipath_analysis")



def read_pickle(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data



def test_GNSS_MultipathAnalysis_sp3_file():
    """
    Test the results from OPEC2022 and sp3 files
    """
    rinObs_file =  r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\NMBUS_SAMSUNG_S20.20o"
    sp3Nav_file = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\SP3\NMBUS_2020 10 30.SP3"
    expected_res = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\Multipath_analysis\Tests\analysisResults_NMBUS.pkl"
    result = GNSS_MultipathAnalysis(rinObsFilename=rinObs_file, sp3NavFilename_1=sp3Nav_file,
                                    plotEstimates=False,
                                    plot_polarplot=False)
    
    expected_result = read_pickle(expected_res)
    # Compare the result with the expected result
    assert_almost_equal(expected_result['Sat_position']['G']['Position'][2],  result['Sat_position']['G']['Position'][2],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['G']['Position'][4],  result['Sat_position']['G']['Position'][4],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['G']['Position'][30],  result['Sat_position']['G']['Position'][30],decimal=4)

    assert_almost_equal(expected_result['Sat_position']['G']['Azimut'],  result['Sat_position']['G']['Azimut'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['Azimut'],  result['Sat_position']['E']['Azimut'],decimal=4)

    assert_almost_equal(expected_result['Sat_position']['G']['Elevation'],  result['Sat_position']['G']['Elevation'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['Elevation'],  result['Sat_position']['E']['Elevation'],decimal=4)


def test_GNSS_MultipathAnalysis_broadcast_navfile():
    """
    Test the results from OPEC2022 and use of broadcast ephemerids
    """
    rinObs_file =  r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\ObservationFiles\OPEC00NOR_S_20220010000_01D_30S_MO_3.04.rnx"
    broadNav_file = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\TestData\NavigationFiles\BRDC00IGS_R_20220010000_01D_MN.rnx"
    expected_res = r"C:\Users\perhe\OneDrive\Documents\Python_skript\GNSS_repo\Multipath_analysis\Tests\analysisResults_OPEC00NOR_S_20220010000_01D_30S_MO.pkl"
    result = GNSS_MultipathAnalysis(rinObsFilename=rinObs_file, broadcastNav1=broadNav_file,
                                    plotEstimates=False,
                                    plot_polarplot=False)
    expected_result = read_pickle(expected_res)
    # Compare the result with the expected result
    assert_almost_equal(expected_result['Sat_position']['G']['Position'][1],  result['Sat_position']['G']['Position'][1],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['G']['Position'][8],  result['Sat_position']['G']['Position'][8],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['G']['Position'][10], result['Sat_position']['G']['Position'][10],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['R']['Position'][1],  result['Sat_position']['R']['Position'][1],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['R']['Position'][7],  result['Sat_position']['R']['Position'][7],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['R']['Position'][8],  result['Sat_position']['R']['Position'][8],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['Position'][1],  result['Sat_position']['E']['Position'][1],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['Position'][3],  result['Sat_position']['E']['Position'][3],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['Position'][8],  result['Sat_position']['E']['Position'][8],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['C']['Position'][5],  result['Sat_position']['C']['Position'][5],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['C']['Position'][6],  result['Sat_position']['C']['Position'][6],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['C']['Position'][9],  result['Sat_position']['C']['Position'][9],decimal=4)

    assert_almost_equal(expected_result['Sat_position']['G']['Azimut'],  result['Sat_position']['G']['Azimut'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['R']['Azimut'],  result['Sat_position']['R']['Azimut'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['Azimut'],  result['Sat_position']['E']['Azimut'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['C']['Azimut'],  result['Sat_position']['C']['Azimut'],decimal=4)

    assert_almost_equal(expected_result['Sat_position']['G']['Elevation'],  result['Sat_position']['G']['Elevation'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['R']['Elevation'],  result['Sat_position']['R']['Elevation'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['E']['Elevation'],  result['Sat_position']['E']['Elevation'],decimal=4)
    assert_almost_equal(expected_result['Sat_position']['C']['Elevation'],  result['Sat_position']['C']['Elevation'],decimal=4)

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

# Run the tests
if __name__ == '__main__':
    pytest.main()


    
