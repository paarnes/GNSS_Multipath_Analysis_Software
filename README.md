# GNSS Multipath Analysis
[![Python application](https://github.com/paarnes/GNSS_Multipath_Analysis_Software/actions/workflows/run-tests.yml/badge.svg)](https://github.com/paarnes/GNSS_Multipath_Analysis_Software/actions/workflows/run-tests.yml)

GNSS_MultipathAnalysis is a software for analyzing the multipath effect on Global Navigation Satellite Systems (GNSS). This software is largely based on the MATLAB software [GNSS_Receiver_QC_2020](https://gitlab.com/bjro/GNSS_reading_protocol/-/tree/main/GNSS_Receiver_QC_2020) made by Bjørn-Eirik Roald. Mainly it follows the same logic, just with Python syntax. However, there are added some other features like for instance:
* Possible to use broadcasted ephemerides (not only SP3 files)
* Also support RINEX v2.xx observation files 
* Makes polar plot of each satellite for each system
* Makes polar plot that shows the multipath effect as function of azimuth and elevation angle.
* Plots the Signal-To-Noise Ratio (SNR) wrt to time and elevation angle
* Makes polar plot that shows the Signal-To-Noise Ratio (SNR) as function of azimuth and elevation angle
* Possible to choose which navigation system to run analysis on (not hardcoded anymore)
* Summary of the number of cycle slips detected in total (both ionospheric residuals and code- phase difference)

A considerable part of the results has been validated by comparing the results with estimates from RTKLIB. This software will be further developed, and feedback and suggestions are therefore gratefully received. Don't hesitate to report if you find bugs or missing functionality. Either by e-mail or by raising an issue here in GitHub.

The main function is called "GNSS_MultipathAnalysis.py" and takes in different arguments. Two arguments are mandatory:
* A RINEX Observation file
* A sp3/eph file containing precise satellite coordinates or a RINEX 3 navigation file 

The rest of the arguments are optional. Their default values are described in the function description. By default, this software will provide the results in forms of plots and an analysis report as a text file. In addition, it exports the results as a pickle file which can be imported as a dictionary in python for instance.

## Installation 
To install the required packages, run:     
`pip install -r requirements.txt`   
where the *requirements.txt* is located [here](https://github.com/paarnes/GNSS/blob/master/src/requirements.txt).

Note: In the example plots, TEX is used to get prettier text formatting. However, this requires TEX/LaTex to be installed on your computer. The program will first try to use TEX, and if it's not possible, standard text formatting will be used. So TEX/LaTex is not required to run the program and make plots.

## The steps are:
1. Reads in the RINEX observation file 
2. Reads the RINEX navigation file or the precise satellite coordinates in SP3-format (depends on what’s provided)
3. If a navigation file is provided, the satellite coordinates will be transformed from Kepler-elements to ECEF for GPS, Galileo and BeiDou. For GLONASS the navigation file is containing a state vector. The coordinates then get interpolated to the current epoch by solving the differential equation using a 4th order Runge-Kutta. If a SP3 file is provided, the interpolation is done by a barycentric Lagrange interpolation. 
4. Satellites elevation and azimuth angles get computed. 
5. Cycle slip detection by using both ionospheric residuals and a code-phase combination. These linear combinations are given as
$$\dot{I} = \frac{1}{\alpha-1}\left(\Phi_1 - \Phi_2\right)/\Delta t$$
$$d\Phi_1R_1 = \Phi_1 - R_1$$
 The threshold values can be set by the user, and the default values are set to $0.0667 [\frac{m}{s}]$ and $6.67[\frac{m}{s}]$ for the ionospheric residuals and code-phase combination respectively. 

6. Multipath estimates get computed by making a linear combination of the code and phase observation. PS: A dual frequency receiver is necessary because observations from two different bands/frequency are needed. 
$$MP_1 = R_1 - \left(1+\frac{2}{\alpha - 1}\right)\Phi_1 + \left(\frac{2}{\alpha - 1}\right)\Phi_2$$
where $R_1$ is the code observation on band 1, $\Phi_1$ and $\Phi_2$ is phase observation on band 1 and band 2 respectively. Furthermore $\alpha$ is the ratio between the two frequency squared $\alpha=\frac{{f}^2_1}{{f}^2_2}$
7. Based on the multipath estimates computed in step 6, both weighted and unweighted RMS-values get computed. The RMS value is given as
$$RMS=\sqrt{\frac{\sum\limits_{i=1}^{N_{sat}}\sum\limits_{j=1}^{N_{epohcs}} MP_{ij}}{N_{est}}}$$ 
For the weighted RMS value, the satellite elevation angle is used in a weighting function defined as $$w =\frac{1}{4sin^2\beta}$$ for every estimates with elevation angle $\beta$ is below $30^{\circ}$ and $w =1$ for $\beta > 30^{\circ}$. 
8. Several plot will be generated (if not set to FALSE): 
    * Ionospheric delay wrt time and zenith mapped ionospheric delay (combined)
		<p align="center">
			<img src="https://github.com/paarnes/GNSS/blob/master/Results_example/2022_30sec/Graphs/Galileo_ionospheric_delay_combined.png" width="630"/>
		</p>	
    * The Multipath effect plotted wrt time and elevation angle (combined)
		<p align="center">
			<img src="https://github.com/paarnes/GNSS/blob/master/Results_example/2018_1sec/Graphs/GLONASS_C1C_C2P_MP_combined.png" width="630"/>
		</p>	
    * Barplot showing RMS values for each signal and system
		<p align="center">
			<img src="https://github.com/paarnes/GNSS/blob/master/Results_example/2022_30sec/Graphs/Barplot_RMS_all.png" width="630"/>
		</p>
    * Polar plot of the multipath effect as function of elevation angle and azimuth
		<p align="center">
			<img src="https://github.com/paarnes/GNSS/blob/master/Results_example/2022_30sec/Graphs/MP_GPS_C1C.png" width="630"/>
		</p>
    * Polar plot of each observed satellite in the system
		<p align="center">
			<img src="https://github.com/paarnes/GNSS/blob/master/Results_example/2022_30sec/Graphs/Skyplot_GLONASS.png" width="630"/>
		</p>

    * Signal-To-Noise Ratio (SNR) plotted wrt time and elevation angle (combine) 
		<p align="center">
			<img src="https://github.com/paarnes/GNSS/blob/master/Results_example/2022_30sec/Graphs/SNR_GPS_S1C.png" width="630"/>
		</p>

    * Polar plot of Signal-To-Noise Ratio (SNR)
		<p align="center">
			<img src="https://github.com/paarnes/GNSS/blob/master/Results_example/2018_1sec/Graphs/SNR_Polar_GPS_S2W.png" width="630"/>
		</p>
9. Exporting the results as a pickle file which easily can be imported into python as a dictionary
10. The results in form of a report get written to a text file with the same name as the RINEX observation file. 


## How to run it
An example file on how to call the program is located [here](https://github.com/paarnes/GNSS_Multipath_Analysis_Software/blob/master/src/Examples_on_how_to_run_it.ipynb). This will show some examples
on how to run the analysis with different user defined arguments, how to read in the resultfile (pickle file), and in addition it shows how to use only the RINEX
reading routine. The most simple example on how to run the code is:

```
from GNSS_MultipathAnalysis import GNSS_MultipathAnalysis  

rinObs_file = 'OPEC00NOR_S_20220010000_01D_30S_MO_3.04' 
SP3_file    = 'SP3_20220010000.eph'   
analysisResults = GNSS_MultipathAnalysis(rinex_obs_file, sp3NavFilename_1=SP3_file)
```


