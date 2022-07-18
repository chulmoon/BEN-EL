# Bayesian Elastic Net based on Empirical Likelihood

We propose a Bayesian elastic net that uses empirical likelihood and develop an efficient tuning of Hamiltonian Monte Carlo for posterior sampling. The proposed model relaxes the assumptions on the identity of the error distribution, performs well when the variables are highly correlated, and enables more straightforward inference by providing posterior distributions of the regression coefficients. Simulation studies and case study on air pollution data are carried and show that the proposed method performs better than the existing methods.

# Simulation Study
 - Generating data
   * Run R files in `./simulation/data/` to generate data for simulation studies.
 - Fitting models
   * For BEN-EL, the initial step size ($\epsilon$) and penalty parameters ($\lambda_1$ and $\lambda_2$) are estimated first and fit the model. For example, for simulation 1, run `simulation1_BENEL_parameter.R` first and fit the model using `simulation1_BENEL.R`.
   * For the other methods (BEN, BL, EN, and LADL), run the corresponding R files. For example, for the BEN model of simulation 1, run `simulation1_BEN.R`.
   * For summary of results, run `simulation1_result.R`, `simulation2_result.R`, and `simulation3_result.R`.

# Air Pollution Case Study

- Run `pollution.R` for applications and plots. `pollution.Rdata` is the data file from McDonald et al. [^1].

[^1]: McDonald, G. C. and Schwing, R. C. (1973) Instabilities of regression estimates relating air pollution to mortality, Technometrics, 15, 463-482.
