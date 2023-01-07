# GNSS_MultipathAnalysis
GNSS_MultipathAnalysis is a software for analyzing the multipath effect on Global Navigation Satellite Systems (GNSS). This software is largely based on the MATLAB software "GNSS_Receiver_QC_2020" made by Bjørn Eirik Roald. Mainly it follows the exact same logic, just with Python syntax instead. Some other features is added like:
* Possible to use broadcasted ephemerides (not only SP3 files)
* Makes skyplot that shows the multipath effect as function of azimut and elevation angle.
* Possible to choose which navigation system to run analysis on (no hardcoded anymore)
