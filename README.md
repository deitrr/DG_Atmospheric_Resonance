[![DOI](https://zenodo.org/badge/793829283.svg)](https://zenodo.org/doi/10.5281/zenodo.11154642)

# DG_Atmospheric_Resonance
Scripts for plots and analysis of Deitrick &amp; Goldblatt atmospheric resonance paper: [To be updated]

The simulation output necessary to plot all figures is located here: [https://doi.org/10.20383/103.0960](https://doi.org/10.20383/103.0960)

## Requirements

This code requires Python 3 as well as the following Python libraries: numpy, scipy, matplotlib, netCDF4, cdo, basemap, pyshtools

## Post-processing and plotting

Outputs from the atmosphere model were monthly history files ("h0") for 60 years of integration and daily history files ("h1") with hourly data for an additional year, usually year 61. These history files are not included in the data repository above as the final size would be far too large to share. 

The original history files are processed using CDO into merged and averaged netCDF4 files that contain all fields or only some that I needed for analysis. In the data repository, the merged files, which are very large, are retained for the simulations used in the main analysis in the directory `sims_main`. For other simulations, only the averaged (and much smaller) netCDF4 files are retained in the directory `sims_suppl`. 

The merged and averaged files are processed using `merge_h0.py`, `merge_h1.py`, and `time_mean_files.py`. These are executed for a single simulation in the example PBS files `proc_merge_h0.job` and `proc_merge_h1.job`, which can be turned into simple bash scripts by removing the PBS preamble. These files are included here for record keeping, i.e., to transparently document the process I used. Since the original history files are not included in the data repository, these scripts are not needed for any plotting or post-processing.

The most intensive part of post-processing is what I called "diurnal-averaging around the Solar longitude". This is done using `diurnal_solar_avg.py`, which can be simply executed from the command line, though I don't recommend running it as is because it will take a very long time and will use up to 20 GB of memory. It is better to open that file, comment out certian items in the "fields" list, and run in piecemeal, or else to use a scheduler on a cluster. The output from that script is also included in the data repo, so it is not necessary to run it, but I include it here for transparency and record keeping. It will not process the data for all simulations as is--additional simulations need to be added to the script for that. 

The scripts `figure_x.py` create the plots from the paper, as one would expect. These are included here for reproducibility, transparency, and the odd chance they might be helpful to someone, somewhere. 
