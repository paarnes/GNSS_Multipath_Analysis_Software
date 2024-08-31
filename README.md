# GNSS Multipath Analysis

[![Python application](https://github.com/paarnes/GNSS_Multipath_Analysis_Software/actions/workflows/run-tests.yml/badge.svg)](https://github.com/paarnes/GNSS_Multipath_Analysis_Software/actions/workflows/run-tests.yml)
[![PyPI version](https://badge.fury.io/py/gnssmultipath.svg)](https://badge.fury.io/py/gnssmultipath)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Installing LaTeX (optional)](#installing-latex-optional)
- [Usage](#how-to-run-it)
- [Examples](#some-simple-examples-on-how-to-use-the-software)
- [Compatibility](#compatibility)
- [Documentation](#documentation)
- [License](#license)

## Introduction

GNSS Multipath Analysis is a software tool for analyzing the multipath effect on Global Navigation Satellite Systems (GNSS). It is primarily based on the MATLAB software [GNSS_Receiver_QC_2020](https://gitlab.com/bjro/GNSS_reading_protocol/-/tree/main/GNSS_Receiver_QC_2020) by Bjørn-Eirik Roald, but has been adapted to Python and includes additional features.

## Features
- Support for broadcasted ephemerides (not only SP3 files)
- Compatibility with RINEX v2.xx observation files
- Generates various plots (polar, SNR vs. time/elevation, etc.)
- Extracts GLONASS FCN from RINEX navigation files
- Cycle slip detection and multipath effect estimation
- Results export to CSV and Pickle formats
- Option to choose which navigation system to analyze
- [More...](#the-steps-are)

## Installation

To install the software, use `pip`:

```bash
pip install gnssmultipath
```

### Prerequisites
- **Python >=3.8**: Ensure you have Python 3.8 or newer installed.
- **LaTeX** (optional): Required for generating plots with LaTeX formatting.

### Installing LaTeX (optional)
- On Ubuntu: `sudo apt-get install texlive-full`
- On Windows: Download and install from [MiKTeX](https://miktex.org/download)
- On MacOS: `brew install mactex`

## How to Run It

To run the GNSS Multipath Analysis, follow these steps:

1. **Prepare your data:**
   - Ensure you have RINEX observation and navigation files.
   
2. **Run the analysis:**
   - Use the following code snippet in your Python script or notebook:

```python
from gnssmultipath import GNSS_MultipathAnalysis

rinObs_file = 'your_observation_file.XXO'
rinNav_file = 'your_navigation_file.XXN'
results = GNSS_MultipathAnalysis(rinex_obs_file=rinObs_file, sp3NavFilename_1=rinNav_file)
```

3. **Analyze the results:**
   - Results will be generated in the specified output directory in the form of plots and text reports.

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

where $R_1$ is the code observation on band 1, $\Phi_1$ and $\Phi_2$ is phase observation on band 1 and band 2 respectively. Furthermore $\alpha$ is the ratio between the two frequency squared 

$$\alpha=\frac{{f}^2_1}{{f}^2_2}$$

7. Based on the multipath estimates computed in step 6, both weighted and unweighted RMS-values get computed. The RMS value has unit _meter_, and is given by

$$RMS=\sqrt{\frac{\sum\limits_{i=1}^{N_{sat}}\sum\limits_{j=1}^{N_{epohcs}} MP_{ij}}{N_{est}}}$$

For the weighted RMS value, the satellite elevation angle is used in a weighting function defined as 

$$w =\frac{1}{4sin^2\beta}$$ 

for every estimates with elevation angle $\beta$ is below $30^{\circ}$ and $w =1$ for $\beta > 30^{\circ}$.

8. Several plot will be generated (if not set to FALSE):
    * Ionospheric delay wrt time and zenith mapped ionospheric delay (combined)
		<p align="center">
			<img src="https://github.com/paarnes/GNSS_Multipath_Analysis_Software/blob/master/Results_example/2022_30sec/Graphs/Galileo_ionospheric_delay_combined.png?raw=true" width="630"/>
		</p>
    * The Multipath effect plotted wrt time and elevation angle (combined)
		<p align="center">
			<img src="https://github.com/paarnes/GNSS_Multipath_Analysis_Software/blob/master/Results_example/2018_1sec/Graphs/GLONASS_C1C_C2P_MP_combined.png?raw=true" width="630"/>
		</p>
    * Barplot showing RMS values for each signal and system
		<p align="center">
			<img src="https://github.com/paarnes/GNSS_Multipath_Analysis_Software/blob/master/Results_example/2022_30sec/Graphs/Barplot_RMS_all.png?raw=true" width="630"/>
		</p>
    * Polar plot of the multipath effect as function of elevation angle and azimuth
		<p align="center">
			<img src="https://github.com/paarnes/GNSS_Multipath_Analysis_Software/blob/master/Results_example/2022_30sec/Graphs/MP_GPS_C1C.png?raw=true" width="630"/>
		</p>
    * Polar plot of each observed satellite in the system
		<p align="center">
			<img src="https://github.com/paarnes/GNSS_Multipath_Analysis_Software/blob/master/Results_example/2022_30sec/Graphs/Skyplot_GLONASS.png?raw=true" width="630"/>
		</p>

    * Signal-To-Noise Ratio (SNR) plotted wrt time and elevation angle (combine)
		<p align="center">
			<img src="https://github.com/paarnes/GNSS_Multipath_Analysis_Software/blob/master/Results_example/2022_30sec/Graphs/SNR_GPS_S1C.png?raw=true" width="630"/>
		</p>

    * Polar plot of Signal-To-Noise Ratio (SNR)
		<p align="center">
			<img src="https://github.com/paarnes/GNSS_Multipath_Analysis_Software/blob/master/Results_example/2018_1sec/Graphs/SNR_Polar_GPS_S2W.png?raw=true" width="630"/>
		</p>
9. Exporting the results as a pickle file which easily can be imported into python as a dictionary
10. The results in form of a report get written to a text file with the same name as the RINEX observation file.
11. The estimated values are also written to a CSV file by default
		<p align="center">
			<img src="https://github.com/paarnes/GNSS_Multipath_Analysis_Software/blob/master/Results_example/result_table.png?raw=true" width="630"/>
		</p>

## Some simple examples on how to use the Software:

### Run a multipath analysis using a SP3 file and only mandatory arguments
```python
from gnssmultipath import GNSS_MultipathAnalysis

rinObs_file = 'OPEC00NOR_S_20220010000_01D_30S_MO_3.04'
SP3_file    = 'SP3_20220010000.eph'
analysisResults = GNSS_MultipathAnalysis(rinex_obs_file=rinObs_file, sp3NavFilename_1=SP3_file)
```

### Run a multipath analysis using a RINEX navigation file with SNR and a defined datarate for ephemerides
```python
from gnssmultipath import GNSS_MultipathAnalysis

# Input arguments
rinObs_file = 'OPEC00NOR_S_20220010000_01D_30S_MO_3.04'
rinNav_file = 'BRDC00IGS_R_20220010000_01D_MN.rnx'
output_folder = 'C:/Users/xxxx/Results_Multipath'
cutoff_elevation_angle = 10  # drop satellites lower than 10 degrees
nav_data_rate = 60  # desired datarate for ephemerides (to improve speed)

analysisResults = GNSS_MultipathAnalysis(rinex_obs_file=rinObs_file,
                                         broadcastNav1=rinNav_file,
                                         include_SNR=True,
                                         outputDir=output_folder,
                                         nav_data_rate=nav_data_rate,
                                         cutoff_elevation_angle=cutoff_elevation_angle)
```

### Read a RINEX observation file
```python
from gnssmultipath import readRinexObs

rinObs_file = 'OPEC00NOR_S_20220010000_01D_30S_MO_3.04'
GNSS_obs, GNSS_LLI, GNSS_SS, GNSS_SVs, time_epochs, nepochs, GNSSsystems, \
        obsCodes, approxPosition, max_sat, tInterval, markerName, rinexVersion, recType, timeSystem, leapSec, gnssType, \
        rinexProgr, rinexDate, antDelta, tFirstObs, tLastObs, clockOffsetsON, GLO_Slot2ChannelMap, success = \
        readRinexObs(rinObs_file)
```

### Read a RINEX navigation file (v.3)
```python
from gnssmultipath import Rinex_v3_Reader

rinNav_file = 'BRDC00IGS_R_20220010000_01D_MN.rnx'
navdata = Rinex_v3_Reader().read_rinex_nav(rinNav_file, data_rate=60)
```

### Read in the results from an uncompressed Pickle file
```python
from gnssmultipath import PickleHandler

path_to_picklefile = 'analysisResults.pkl'
result_dict = PickleHandler.read_pickle(path_to_picklefile)
```

### Read in the results from a compressed Pickle file
```python
from gnssmultipath import PickleHandler

path_to_picklefile = 'analysisResults.pkl'
result_dict = PickleHandler.read_zstd_pickle(path_to_picklefile)
```

## Compatibility
- **Python Versions:** Compatible with Python 3.8 and above.
- **Dependencies:** All dependencies will be automatically installed with `pip install gnssmultipath`.

## Documentation

Detailed documentation is available [here](https://gnssmultipath.readthedocs.io/).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
