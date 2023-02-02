# GNSS Multipath Analysis
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

## The steps are:
1. Reads in the RINEX observation file 
2. Reads the RINEX navigation file or the precise satellite coordinates in SP3-format (depends on whats provided)
3. If a navigation file is provided, the satellite coordinates will be transformed from Kepler-elemements to ECEF for GPS,Galileo and BeiDou. For GLONASS the navigation file is containing a state vector. The coordinates then get interpolated to the current epoch by solving the differential equation using a 4th order Runge-Kutta. If a SP3 file is provided, the interpolation is done by a barycentric lagrange interpolation. 
4. Satellites elevation and azimut angles get computed. 
5. Cycle slip detection by using both ionospheric residuals and a code-phase combination. These linear combination is given as
$$\dot{I} = \frac{1}{\alpha-1}\left(\Phi_1 - \Phi_2\right)/\Delta t$$
$$d\Phi_1R_1 = \Phi_1 - R_1$$
 The thresholdvalue can be set by user, but the default values are set to $0.0667\: m/s$ and $6.67\: m/s$ for the ionospheric residuals and code-phase combination respectivly. 

6. Multipath estimates get computed by making a linear combination of the code and phase observation. PS: A dual frequency reciver is necessery. Must have observation from two different bands/frequency. 
$$MP_1 = R_1 - \left(1+\frac{2}{\alpha - 1}\right)\Phi_1 + \left(\frac{2}{\alpha - 1}\right)\Phi_2$$
where $R_1$ is the code observation on band 1, $\Phi_1$ and $\Phi_2$ is phase observation on band 1 and band 2 repectivly. Furthermore $\alpha$ is the ratio between the two frequenzy squered $\alpha=\frac{{f}^2_1}{{f}^2_2}$
7.Based on the multipath estimates computed in step 6, both weighted and unweighted RMS-values get computed. The RMS value is given as
$$RMS=\sqrt{\frac{\sum\limits_{i=0}^{N_{sat}}\sum\limits_{j=0}^{N_{epohcs}} MP_{ij}}{N_{est}}}$$ 
For the weighted RMS value, the satellite elevation angle is used in a weighting function $$w =\frac{1}{4sin^2\beta}$$ for every estimates with elevation angle $\beta$ is below $30^{\circ}$. Else $w =1$. 
8. Several plot will be generated (if not set to FALSE): 
    * Ionospheric delay wrt time and zenit mapped ionospheric delay (combined)
    * THe Multipath effect plotted wrt time and elevation angle (combined)
    * Barplot showing RMS values for each signal and system
    * Polar plot of the multipath effect as funciton of elevaion angle and azimut
    * Polar plot of each observed satellite in the system
9. Exporting the results as a pickle file which easily can be imported into python as a dictionary
10. The results in form of a report get written to a text file with the same name as the RINEX observation file. 




