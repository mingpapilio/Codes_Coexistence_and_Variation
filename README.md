# Codes_Coexistence_and_Variation
This repository is for open access to source codes and data for the manuscript "Antagonistic effects of long- and short-term 
environmental variation on species coexistence".

1. **PLEASE READ** the notes at the beginning of each files, containing the important descriptions of key parameters in simulations and notes for plotting. 

2. In "codes" folder, "pc2_tseries.c" generates time series data while "rpc2_avrt.c" generates the data for the heat map of coexistence. Time series data is used in Fig. 1, 3, and 4. Heat map data is used in Fig. 2 and 5. 

3. In "data" folder, raw data are sorted according to figure and panel indices. Some folders have R files to help plotting. 

4. In "codes" folder, please **DO NOT** change the relative location of these files because of the random number generators files (dSFMT-src-2.2.3). The environments for simulations were (1) clang-800.0.42.1 (Mac OSX 10.12.3), (2) gcc version 5.4.0 20160609 (Ubuntu 5.4.0-6ubuntu1~16.04.10), and (3) gcc version 9.2.0 (Homebrew on Mac OSX). 

Please email ming.liu@zoo.ox.ac.uk if there is any problem!
(Ming Liu)
