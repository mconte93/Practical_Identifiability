# Practical_Identifiability_ConteEtAl2024
This repository contains all the data and codes used to perform simulations as described in

Conte, Martina, et al. "Structural and practical identifiability of contrast transport models for DCE-MRI." PLOS Computational Biology 20.5 (2024): e1012106.

The published version of the paper can be found at the link [https://doi.org/10.1101/2023.12.19.572294](https://doi.org/10.1371/journal.pcbi.1012106).
For details we refer the interested reader to this publication. If you use this software in your work then please cite the above named paper.


### Copyright notice

Practical_Identifiability_ConteEtAl2024: Matlab code to perform the analysis on practical identifiability on contrast transport models for DCE-MRI. 
Copyright (C) 2024 M. Conte

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.

## Organization of the repository
This repository is organized as follows.

 - "Data" folder contains the original scans used for the GBM and Breast cancer cases.
   
 - "MLE_LTK" folder contains the codes used for the results obtained with the LTK model. It is divided into 6 subfolders: AA study, RA study, RR GBM study, RR Breast study, noise study, and smoothing study. Each of them contains data and codes to replicate the results shown in the paper.

 - "MLE_PM", "MLE_TK", "MLE_ETK", and "MLE_mLTK" folders contain the codes used for the results obtained with the PM, TK, ETK, and mLTK models, respectively. Each of them is divided into 4 subfolders for AA study, RA study, RR GBM study, and RR Breast study, and they contain data and codes to replicate the results shown in the supplementary material of the paper.
