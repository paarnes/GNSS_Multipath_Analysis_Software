# GNSS_MultipathAnalysis
GNSS_MultipathAnalysis is a software for analyzing the multipath effect on Global Navigation Satellite Systems (GNSS). This software is largely based on the MATLAB software "GNSS_Receiver_QC_2020" made by Bjørn Eirik Roald. Mainly it follows the same logic, just with Python syntax. However there are added some other features like for instance:
* Possible to use broadcasted ephemerides (not only SP3 files)
* Makes polar plot of each satellite for each system
* Makes polar plot that shows the multipath effect as function of azimut and elevation angle.
* Possible to choose which navigation system to run analysis on (no hardcoded anymore)
* Summary of the number of cycle slips detected i total (both ionospheric residuals and code- phase difference)

The main function is called "GNSS_MultipathAnalysis.py" and takes in different arguments. Two arguments are mandatory:
* A RINEX 3 Observation file
* A sp3/eph file containing precise satellite coordinates or a RINEX 3 navigation file 

The rest of the arguments are optional. Their default values are described in the function description. By default this software will provide the results in forms of plots and a analysis report as a text file. In addition it exports the results as a pickle file which can be imported as a dictionary in python for instance. 

The steps are:
1. Reads in the RINEX observation file 
2. Reads the RINEX navigation file or the precise satellite coordinates in SP3-format (depends on whats provided)
3. If a navigation file is provided, the satellite coordinates will be transformed from Kepler-elemements to ECEF for GPS,Galileo and BeiDou. For GLONASS the navigation file is containing a state vector. The coordinates then get interpolated to the current epoch by solving the differential equation using a 4th order Runge-Kutta. If a SP3 file is provided, the interpolation is done by a barycentric lagrange interpolation. 
4. Satellites elevation and azimut angles get computed. 
5. Cycle slip detection by using both ionospheric residuals and a code-phase combination. (LEGGE TIL FORMLER?). The thresholdvalue can be set by user.
6. Multipath estimates get computed, both weighted and unweighted RMS-values. For the weighted RMS value, the satellite elevation angle is used in a weighting function (LEGGE TIL HER?)


