# Codes_Coexistence_and_Variation
This repository is for open access to source codes and data for the manuscript "Antagonistic effects of long- and short-term 
environmental variation on species coexistence".

1. **PLEASE READ** the notes at the beginning of each files, containing the important descriptions of key parameters in simulations and notes for plotting. 

2. In **/Codes** folder, there are two subfolders for running and a random number generator for C (dSFMT-src-2.2.3). The subfolder 'TimeSeries' contains the codes for generating all time series data used in the paper; while the subfolder 'HeatMap' contains the codes for generating all of the coexistence heat map results. Note that the relative location of folders are important, including the dSFMT folder, and gen_beta files are the normal-distributed random number generator. Lastly, sample running command and key variables are explained at the beginning of the main code file (pc2_tseries.c and rpc2_vart.c), please check them before running. 

3. In **/Data** folder, raw data are sorted according to figure and panel indices. Some folders have R files to help plotting. 

Please email ming.liu@zoo.ox.ac.uk if you have any question! (Ming Liu)
