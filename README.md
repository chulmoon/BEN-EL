# Bayesian Elastic Net based on Empirical Likelihood

We propose a Bayesian elastic net that uses empirical likelihood and develop an efficient tuning of Hamiltonian Monte Carlo for posterior sampling. The proposed model relaxes the assumptions on the identity of the error distribution, performs well when the variables are highly correlated, and enables more straightforward inference by providing posterior distributions of the regression coefficients. Simulation studies and case study on air pollution data are carried and show that the proposed method performs better than the existing methods.

# Data Availability
All data used in simulation studies are generated randomly and the air pollution data are from McDonald et al. [^1]. 

# Simulation and Case Study
1. Simulation studies

 - Generating data
   * Run R files in `./simulation/data/` to generate data for simulation studies.
 - Fitting models
   * For BEN-EL, estimate the parameters first (`simulation1_BENEL_parameter.R`, `simulation2_BENEL_parameter.R`, `simulation3_BENEL_l1_parameter.R` , and `simulation3_BENEL_l2_parameter.R`) and fit the models (`simulation1_BENELr.R`, `simulation2_BENELr.R`, `simulation3_BENEL_l1.R` , and `simulation3_BENEL_l2.R`)
   * For all other models (BEN, BL, EN, and LADL), run the corresponding R files.
   * For summary of results, run `simulation1_result.R`, `simulation2_result.R`, and `simulation3_result.R`.

2. Air pollution case study

- Run `pollution.R` for applications and plots. `pollution.Rdata` is the data file from McDonald et al. [1].

[^1]: McDonald, G. C. and Schwing, R. C. (1973) Instabilities of regression estimates relating air pollution to mortality, Technometrics, 15, 463-482.
