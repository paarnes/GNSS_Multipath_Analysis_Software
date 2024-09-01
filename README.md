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
- [How to run it](#how-to-run-it)
- [Compatibility](#compatibility)
- [License](#license)
- [User arguments](#user-arguments)
- [Examples](#some-simple-examples-on-how-to-use-the-software)


## Introduction

GNSS Multipath Analysis is a software tool for analyzing the multipath effect on Global Navigation Satellite Systems (GNSS). It is primarily based on the MATLAB software [GNSS_Receiver_QC_2020](https://gitlab.com/bjro/GNSS_reading_protocol/-/tree/main/GNSS_Receiver_QC_2020) by Bjørn-Eirik Roald, but has been adapted to Python and includes additional features.

## Features
- Supports broadcasted ephemerides (RINEX nagiation files) and SP3 files
- Supports both RINEX v2.xx and v3.xx observation files
- Generates various plots like:
	- Ionospheric delay wrt time and zenith mapped ionospheric delay (combined)
	- The Multipath effect plotted wrt time and elevation angle (combined)
	- Barplot showing RMS values for each signal and system
	- Polar plot of the multipath effect and the Signal to nosie ratio (SNR)
	- polarplots of SNR and multipath
	- SNR vs. time/elevation
	-
- Extracts GLONASS FCN from RINEX navigation files
- Cycle slip detection and multipath effect estimation
- Results export to CSV and Pickle formats
- Option to choose which navigation system to analyze


## Installation

To install the software to your Python environment using ``pip``:

```bash
pip install gnssmultipath
```

### Prerequisites
- **Python >=3.8**: Ensure you have Python 3.8 or newer installed.
- **LaTeX** (optional): Required for generating plots with LaTeX formatting.

### Installing LaTeX (optional)
- On Ubuntu: `sudo apt-get install texlive-full`
- On Windows: Download and install from [MiKTeX](https://miktex.org/download)
- On MacOS: ` brew install --cask mactex`

## How to Run It

To run the GNSS Multipath Analysis, import the main function and specify the RINEX observation and navigation/SP3 files you want to use. To perform the analysis with default settings and by using a navigation file:

```python
from gnssmultipath import GNSS_MultipathAnalysis

outputdir = 'path_to_your_output_dir'
rinObs_file = 'your_observation_file.XXO'
rinNav_file = 'your_navigation_file.XXN'
analysisResults = GNSS_MultipathAnalysis(rinObs_file,
                                         broadcastNav1=rinNav_file,
                                         outputDir=outputdir)
```

## The steps are:
1. Reads in the RINEX observation file
2. Reads the RINEX navigation file or the precise satellite coordinates in SP3-format (depends on what’s provided)
3. If a navigation file is provided, the satellite coordinates will be transformed from Kepler-elements to ECEF for GPS, Galileo and BeiDou. For GLONASS the navigation file is containing a state vector. The coordinates then get interpolated to the current epoch by solving the differential equation using a 4th order Runge-Kutta. If a SP3 file is provided, the interpolation is done by a barycentric Lagrange interpolation.
4. Satellites elevation and azimuth angles get computed.
5. Cycle slip detection by using both ionospheric residuals and a code-phase combination. These linear combinations are given as

$$
\dot{I} = \frac{1}{\alpha-1}\left(\Phi_1 - \Phi_2\right)/\Delta t
$$

$$
d\Phi_1R_1 = \Phi_1 - R_1
$$

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
			<img src="https://github.com/paarnes/GNSS_Multipath_Analysis_Software/blob/master/Results_example/result_table.PNG?raw=true" width="830"/>
		</p>

## User arguments

The `GNSS_MultipathAnalysis` function accepts several keyword arguments that allow for detailed customization of the analysis process. Below is a list of the first five arguments:

- **rinObsFilename** (`str`):  
  Path to the RINEX 3 observation file. This is a required argument.

- **broadcastNav1** (`Union[str, None]`, optional):  
  Path to the first RINEX navigation file. Default is `None`.

- **broadcastNav2** (`Union[str, None]`, optional):  
  Path to the second RINEX navigation file (if available). Default is `None`.

- **broadcastNav3** (`Union[str, None]`, optional):  
  Path to the third RINEX navigation file (if available). Default is `None`.

- **broadcastNav4** (`Union[str, None]`, optional):  
  Path to the fourth RINEX navigation file (if available). Default is `None`.

<details>
  <summary>More...</summary>

  - **sp3NavFilename_1** (`Union[str, None]`, optional):  
    Path to the first SP3 navigation file. Default is `None`.

  - **sp3NavFilename_2** (`Union[str, None]`, optional):  
    Path to the second SP3 navigation file (optional). Default is `None`.

  - **sp3NavFilename_3** (`Union[str, None]`, optional):  
    Path to the third SP3 navigation file (optional). Default is `None`.

  - **desiredGNSSsystems** (`Union[List[str], None]`, optional):  
    List of GNSS systems to include in the analysis. For example, `['G', 'R']` to include only GPS and GLONASS. Default is all systems (`None`).

  - **phaseCodeLimit** (`Union[float, int, None]`, optional):  
    Critical limit that indicates cycle slip for phase-code combination in m/s. If set to `0`, the default value of `6.667 m/s` will be used. Default is `None`.

  - **ionLimit** (`Union[float, None]`, optional):  
    Critical limit indicating cycle slip for the rate of change of the ionospheric delay in m/s. If set to `0`, the default value of `0.0667 m/s` will be used. Default is `None`.

  - **cutoff_elevation_angle** (`Union[int, None]`, optional):  
    Cutoff angle for satellite elevation in degrees. Estimates with elevation angles below this value will be excluded. Default is `None`.

  - **outputDir** (`Union[str, None]`, optional):  
    Path to the directory where output files should be saved. If not specified, the output will be generated in a sub-directory within the current working directory. Default is `None`.

  - **plotEstimates** (`bool`, optional):  
    Whether to plot the estimates. Default is `True`.

  - **plot_polarplot** (`bool`, optional):  
    Whether to generate polar plots. Default is `True`.

  - **include_SNR** (`bool`, optional):  
    If set to `True`, the Signal-to-Noise Ratio (SNR) from the RINEX observation file will be included in the analysis. Default is `True`.

  - **save_results_as_pickle** (`bool`, optional):  
    If `True`, the results will be saved as a binary pickle file. Default is `True`.

  - **save_results_as_compressed_pickle** (`bool`, optional):  
    If `True`, the results will be saved as a binary compressed pickle file using zstd compression. Default is `False`.

  - **write_results_to_csv** (`bool`, optional):  
    If `True`, a subset of the results will be exported as a CSV file. Default is `True`.

  - **output_csv_delimiter** (`str`, optional):  
    The delimiter to use for the CSV file. Default is a semicolon (`;`).

  - **nav_data_rate** (`int`, optional):  
    The desired data rate for ephemerides in minutes. A higher value speeds up processing but may reduce accuracy. Default is `60` minutes.

  - **includeResultSummary** (`Union[bool, None]`, optional):  
    Whether to include a detailed summary of statistics in the output file, including for individual satellites. Default is `None`.

  - **includeCompactSummary** (`Union[bool, None]`, optional):  
    Whether to include a compact overview of statistics in the output file. Default is `None`.

  - **includeObservationOverview** (`Union[bool, None]`, optional):  
    Whether to include an overview of observation types for each satellite in the output file. Default is `None`.

  - **includeLLIOverview** (`Union[bool, None]`, optional):  
    Whether to include an overview of LLI (Loss of Lock Indicator) data in the output file. Default is `None`.

  - **use_LaTex** (`bool`, optional):  
    If `True`, LaTeX will be used for rendering text in plots, requiring LaTeX to be installed on your system. Default is `True`.

</details>


### Output

- **analysisResults** (`dict`):  
  A dictionary containing the results of the analysis for all GNSS systems.



## Compatibility
- **Python Versions:** Compatible with Python 3.8 and above.
- **Dependencies:** All dependencies will be automatically installed with `pip install gnssmultipath`.


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.




## Some simple examples on how to use the software:

### Run a multipath analysis using a SP3 file and only mandatory arguments
```python
from gnssmultipath import GNSS_MultipathAnalysis

rinObs_file = 'OPEC00NOR_S_20220010000_01D_30S_MO_3.04'
SP3_file    = 'SP3_20220010000.eph'
analysisResults = GNSS_MultipathAnalysis(rinex_obs_file=rinObs_file, sp3NavFilename_1=SP3_file)
```

### Run a multipath analysis using a RINEX navigation file with SNR, a defined datarate for ephemerides and with an elevation angle cut of at 10°
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


### Run analysis without making plots
```python
from gnssmultipath import GNSS_MultipathAnalysis

rinObs_file = 'OPEC00NOR_S_20220010000_01D_30S_MO_3.04'
SP3_file    = 'SP3_20220010000.eph'
analysisResults = GNSS_MultipathAnalysis(rinex_obs_file=rinObs_file, sp3NavFilename_1=SP3_file, plotEstimates=False)
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


