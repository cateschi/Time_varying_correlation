# Estimation of time-varying state correlations in state space models

This repository contains the R and C++ codes that can be used to replicate the results obtained in the paper "Estimation of time-varying state correlations in state space models" by Caterina Schiavoni, Siem Jan Koopman, Franz Palm, Stephan Smeekes and Jan van den Brakel.


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
* ``d`` is the number of state variables for which a diffuse initialisation is used.
* ``gamma_draw`` is a vector containing the Rao-Blackwellised bootstrap filter estimates of the unbounded correlation parameter ($$\gamma_t$$)
* ``hyper_tan``: if ``TRUE``, it uses the hyperbolic function to bound the correlation parameter, otherwise the ``link`` function.
* ``knots`` is the vector of the observations corresponding to the knots.
* ``len`` is the sample size.
* ``nstates`` is the number of state variables in the model.
* ``opti``: if ``TRUE``, it optimizes the function.
* ``outofsample``: if ``TRUE``, it computes the log-likelihood based on the out-of-sample forecast errors.
* ``par`` is the vector of the inital values for the static parameters.
* ``parP10`` is a large number used for the diffuse initialisation.
* ``se`` is the ``len`` times 5 matrix of the standard errors of the GREG estimates.
* ``y`` is the 6 times ``len`` matrix that contains the observed series.


## C++ script *TV_corr.cpp*

This C++ script contains a first set of functions that are needed in order to perform some mathematical operations. The functions that are instead related to the estimation of the time-varying correlation are:
* ``stratifiedResampling_rcpp`` performs startified resampling.
* ``KF_CC_splines_rcpp`` performs the Kalman filter estimation (and log-likelihood evaluation) of the cubic splines model.
* ``ucminf_rcpp_splines`` maximises the the log-likelihood function evaluated in the ``KF_CC_splines_rcpp`` function.
* ``KF_t_rcpp`` computes the prediction step of the Kalman filter.
* ``boot_filter_CC_rcpp`` performs the Rao-blackwellised bootstrap filter estimation of the state vector of the nonlinear model.


In case of queries you can contact the corresponding author at c[dot]schiavoni[at]maastrichtuniversity[dot]nl.
