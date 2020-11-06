# Estimation of time-varying state correlations in state space models

This repository contains the R and C++ codes that can be used to replicate the empirical results obtained in the paper "Estimation of time-varying state correlations in state space models" by Caterina Schiavoni, Siem Jan Koopman, Franz Palm, Stephan Smeekes and Jan van den Brakel.


## R script *TV_corr.R*

This R script should be run in order to replicate the analyses, since it calls the *TV_corr.cpp* script. 

### Functions
The R script contains the following functions:
* ``Weights`` creates the cubic splines weights.
* ``link`` is a link function that bounds its argument between -1 and 1.
* ``KF_CC_known_corr`` performs the Kalman filter estimation (and log-likelihood evaluation) with known values for the time-varying correlation.
* ``KF_CC_splines`` performs the Kalman filter estimation (and log-likelihood evaluation) of the cubic splines model.
* ``KF_CC_const`` performs the Kalman filter estimation (and log-likelihood evaluation) of the model with a time-constant correlation.
* ``IndInf_CC`` performs the estimation of the static parameters of the nonlinear model by indirect inference.

### Inputs
Some inputs are common to many functions. Below you can find the list (in alphabetical order) of all inputs:
* ``beta_hat`` is the vector with the maximum likelihood estimates of the static parameters of the cubic splines model, except for the estimated parameter vector (phi) used in the cubic splines specification of the time-varying correlation.
* ``d`` is the number of state variables for which a diffuse initialisation is used.
* ``eps_sim`` is a matrix with ``nstates-length(states_noerr)+2`` rows and ``len`` times ``sim_h`` columns, and it contains the generated innovations used in the indirect inference estimation.
* ``gamma_draw`` is the vector containing the Rao-Blackwellised bootstrap filter estimates of the unbounded version of the correlation parameter (gamma).
* ``hyper_tan``: if ``TRUE``, it uses the hyperbolic function to bound the correlation parameter, otherwise the ``link`` function.
* ``k`` is the number of knots.
* ``knots`` is the vector of the observations corresponding to the knots.
* ``len`` is the sample size.
* ``nstates`` is the number of state variables in the model.
* ``opti``: if ``TRUE``, it optimizes the function.
* ``outofsample``: if ``TRUE``, it computes the log-likelihood based on the out-of-sample forecast errors.
* ``par`` is the vector of the (inital) values for the static parameters.
* ``parP10`` is a large number used for the diffuse initialisation.
* ``restricted``:  if ``TRUE``, it restricts the model to have a time-constant correlation.
* ``se`` is the ``len`` times 5 matrix of the standard errors of the GREG estimates.
* ``sim_h`` is the number of simulations used in the indirect inference estimation.
* ``states_noerr`` is a vector of indices for the state variables that do not have an error term in the transition equation.
* ``VARgamma`` is the variance of the estimated (by maximum likelihood) parameter vector (phi) used in the cubic splines specification of the time-varying correlation.
* ``W`` is the ``len`` times ``k`` matrix of cubic splines weights (it is the output of the ``Weights`` function).
* ``y`` is the 6 times ``len`` matrix that contains the observed series.

### Outputs
Since also some outputs are common to many functions, below you find another list with all (alphabetically ordered) outputs:
* ``logl`` is the value of the maximised log-likelihood.
* ``opt.fun`` is the value of the minimised objective function used in the indirect inference estimation.
* ``Ptt`` is a list of matrices corresponding to the estimated (co)variances of the Kalman filter estimates of the state variables.
* ``Pttm1`` is a list of matrix corresponding to the predicted (co)variances of the one-step-ahead predictions of the state variables.
* ``st_for`` is the matrix with the standardised forecast errors.
* ``W`` is the ``len`` times ``k`` matrix of cubic splines weights.
* ``xtt`` is the matrix with the Kalman filter estimates of the state variables.
* `xttm1` is the matrix with the one-step-ahead predictions of the state variables.


## C++ script *TV_corr.cpp*

### Functions
This C++ script contains a first set of functions that are needed in order to perform some mathematical operations. The functions that are instead related to the estimation of the time-varying correlation are:
* ``stratifiedResampling_rcpp`` performs startified resampling.
* ``KF_CC_splines_rcpp`` performs the Kalman filter estimation (and log-likelihood evaluation) of the cubic splines model.
* ``ucminf_rcpp_splines`` maximises the the log-likelihood function evaluated in the ``KF_CC_splines_rcpp`` function.
* ``KF_t_rcpp`` computes the prediction step of the Kalman filter.
* ``boot_filter_CC_rcpp`` performs the Rao-blackwellised bootstrap filter estimation of the state vector of the nonlinear model.

Most of the inputs/outputs of the functions included in the C++ script are the same as the ones listed above. We therefore here report only the addtional ones.

### Inputs
*  ``draw_m`` and ``M`` both corresponds the number of particles used for the Rao-Blackwellised bootstrap filter.
* ``init_gamma`` is the value that is used to initialise the Rao-Blackwellised bootstrap filter of gamma (the unbounded version of the correlation parameter).
* ``init_val`` is the vector of the inital values for the static parameters.
* ``Rsel`` is selection matrix of the innovations in the transition equation.
* ``tau_hat`` is the vector with the indirect inference estimates of the static parameters of the nonlinear model.
* ``w`` is a vector of unstandardised weights calculated for the resampling step of the Rao-Blackwellised bootstrap filter's algorithm.

Notice that in the ``KF_t_rcpp`` function, the inputs ``se``, ``xttm1`` and ``y`` are vectors instead of matrices, and ``Pttm1`` is a matrix instead of a list.

### Outputs
* ``att_BF`` is the vector containing the Rao-Blackwellised bootstrap filter estimates of the unbounded version of the correlation parameter (gamma).
* ``CV`` is the vector containing the values for the coefficient of variation, calculated in the Rao-Blackwellised bootstrap filter's algorithm.
* ``ESS`` is the vector containing the values for the effective sample size, calculated in the Rao-Blackwellised bootstrap filter's algorithm.
* ``Ptt_BF`` is the vector containing the estimated variances of the Rao-Blackwellised bootstrap filter estimates of the unbounded version of the correlation parameter (gamma).
* ``resample_set`` is a vector of indices for the particles that have been resampled.
