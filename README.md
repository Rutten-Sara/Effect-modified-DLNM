# Distributed lag non-linear models with effect modification using Laplacian P-splines
Sara Rutten, Thomas Neyens, Elisa Duarte and Christel Faes

## About this repository
This repository contains the Rcodes used to generate the results from the paper "Distributed lag non-linear models with effect modification using Laplacian P-splines".

## Data
Following datasets are necessary to run the code and are included in the data folder:
| Dataset | Description | Downloaded from |
| --- | --- | --- |


## Rcode
The main code covariate_varying_DLNM_Italy.R contains the Rcode of the data application.   
Moreover, the folder **simulations** includes the Rcode that can be used to recover the simulation results. This folder contains several files:
| File | Description |
| --- | --- | 
| Linear_small_area.R | Simulation code for linear effect modification in small areas |
| Linear_large_area.R | Simulation code for linear effect modification in large areas |
| Smooth_small_area.R | Simulation code for smooth effect modification in small areas |
| Smooth_large_area.R | Simulation code for smooth effect modification in large areas |
| Additional Simulations | Similar code to obtain the results for the negative binomial simulations (Supplementary Materials) |

The folder **functions** include the Rcode that is needed to fit the effect modified DLNM model. The folder contains the following files:
| File | Description |
| --- | --- | 
| DLNM_Laplace_covariate.R | Function to fit the DLNM model with effect modification |
| DLNM_Laplace_covariate_poisson.R | Function called by DLNM_Laplace_covariate.R (fit the DLNM model with poisson distribution) |
| DLNM_Laplace_covariate_NB.R | Function called by DLNM_Laplace_covariate.R (fit the DLNM model with negative binomial distribution) |
| predRR_covariate.R | Function to calculate the estimated RR |
| interpret_coefficients.R | Function to calculate the RRR (relative risk ratio) |
| help_functions | Functions called by DLNM_Laplace_covariate.R, predRR_covariate.R and interpret_coefficients.R |
| af_Laplace.R | Calculate attributable fraction |
| BAM_predict.R | Function to calculate the RR using the BAM/GAM package |
