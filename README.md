
# DG_Atmospheric_Resonance
# DG_Pre-Cambrian_Ozone
Scripts for plots and analysis of Deitrick &amp; Goldblatt atmospheric resonance paper: [To be updated]

The simulation output necessary to plot all figures is located here: [To be updated]

## Requirements

This code requires Python 3 as well as the following Python libraries: numpy, scipy, matplotlib, netCDF4, cdo, basemap, pyshtools

## Post-processing

Outputs from the atmosphere model were monthly history files for 60 years of integration and daily history files with hourly data for an additional year, usually year 61. These history files are not included in the data repository above as the final size would be far too large to share. 

The original history files are processed using CDO into merged and averaged netCDF4 files that contain all fields or only some that I needed for analysis. In the data repository, the merged files, which are very large, are retained for the simulations used in the main analysis in the directory `sims_main`. For other simulations, only the averaged (and much smaller) netCDF4 files are retained in the directory `sims_suppl`. 

The merged and averaged files are processed using `merge_h0.py`, `merge_h1.py`, and `time_mean_files.py`. These are executed for a single simulation in the example PBS files `proc_merge_h0.job` and `proc_merge_h1.job`, which can be turned into simple bash scripts by removing the PBS preamble. These files are included here for record keeping, i.e., to transparently document the process I used. Since the original history files are not included in the data repository, these scripts are not needed for any plotting or post-processing.


