# Global Carbon Assimilation System v2 (GCASv2)
GCASv2 is an atmospheric inversion system described in the paper "Regional CO₂ fluxes from 2010 to 2015 inferred from GOSAT XCO₂ retrievals using a new version of the Global Carbon Assimilation System". This repository allows users to generate inversion results for global and regional carbon fluxes.
## Table of Contents
* Overview
* Repository Structure
* System Requirements
  * Hardware Requirements
  * Software Requirements
* Usage Guide
* License
* Contact
# Overview
GCASv2 utilizes Mozart-4 atmospheric transport model and the ensemble root-mean-square filtering algorithm to assimilate satellite XCO2 retrievals to estimate gridded carbon fluxes. For more details about the GCASv2, please refer to the paper "Regional CO2 fluxes from 2010 to 2015 inferred from  GOSAT XCO2 retrievals using a new version of the Global  Carbon Assimilation System". 
# Repo Contents
* src/: Core processing code for GCASv2.
* asmlt/: Contains the compiled executable program.
* input/: A priori flux data and observation files required for the inversion.
* Mozart4/: The atmospheric transport model used by GCASv2.
* output/: Post-processing scripts for the a posteriori flux data.
* Result/: Directory to store inversion results.
* gcas_start: Main script to initiate the GCASv2 processing workflow.
# System Requirements
## Hardware Requirements
When the number of observations and partitions is small and the inversion duration is short, a regular PC can be used to meet the processing requirements. However, in most cases, we may need a computer with higher hardware configuration for inversion to meet the computational needs. The server used in this study has the following specifications:
* RAM: 62G.
* CPU: Intel(R) Xeon(R) CPU E5-2667 v3 @ 3.20GHz.
* GPU: Matrox Electronics Systems Ltd. MGA G200EH.
## Software Requirements
* Operating Systems Tested:
  * Red Hat Enterprise Linux Server release 6.8 (Santiago).
  * Red Hat Enterprise Linux Server release 7.7 (Maipo).
* Programming Language: Fortran (compiled using Intel Fortran Compiler 13.1.0)
# Usage Guide
* 1.Preparation:
  * Ensure you have the required hardware and software environment.
  * Download or clone the repository.
* 2. Input Data:
  * Place the necessary a priori flux data and observation files into the input/ directory.
* 3. Running the Code:
  * Execute the main script gcas_start to initiate the inversion process.
  * The core processing in the src/ directory will run the inversion workflow.
* 4. Post-Processing:
  * The inversion results can be found in the Result/ directory.
# Usage License
This code is provided strictly for academic research purposes. Commercial use, redistribution, or modifications require prior written permission from the project maintainers. Please contact the project maintainers for authorization if you wish to use the code beyond academic research.
# Contact
For questions, collaboration opportunities, or to request commercial use rights, please contact the project maintainers at:
Email: jiangf@nju.edu.cn
