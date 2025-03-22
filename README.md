# GCASv2
## Contents
* Overview
* Repo Contents
* System Requirements
* User Guide
* License
# Overview
This resource is a demo code for Global Carbon Assimilation System version 2 (GCASv2), Which can be used to generate the results based on the atmospheric inversion method in the paper "Regional CO2 fluxes from 2010 to 2015 inferred from  GOSAT XCO2 retrievals using a new version of the Global  Carbon Assimilation System".
# Repo Contents
* src: The folder restores the core processing code of GCASv2. 
* asmlt : The folder restores the compiled executable program.
* input: The folder is used to restore a priori fluxes and observations needed for the inversion.
* Mozart4: Contains the Atmosphere transport model Mozart4.
* output: Contains code that processes the a posteriori flux.
* Result: where is used to restore the inversion results.
* gcas_start: The runing srcipt of GCASv2.
# System Requirements
## Hardware Requirements
When the number of observations and partitions is small and the inversion duration is short, a regular PC can be used to meet the processing requirements. However, in most cases, we may need a computer with higher hardware configuration for inversion to meet the computational needs. The server used in this study has the following specifications:
* RAM: 62G.
* CPU: Intel(R) Xeon(R) CPU E5-2667 v3 @ 3.20GHz.
* GPU: Matrox Electronics Systems Ltd. MGA G200EH.
## Software Requirements
* The code has been tested on the following OS systems:
  * Red Hat Enterprise Linux Server release 6.8 (Santiago).
  * Red Hat Enterprise Linux Server release 7.7 (Maipo).
* The code is implemented in Fortran (Intel Fortran Compiler 13.1.0).
# Usage License
This code is provided for academic research purposes only and may not be used for commercial purposes. Any commercial use, redistribution, or modification of the code requires prior written permission from the project maintainers.
If you wish to use the code, please contact the project maintainers to obtain authorization.

