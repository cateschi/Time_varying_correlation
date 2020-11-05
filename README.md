# Estimation of time-varying state correlations in state space models

This repository contains the R and C++ codes that can be used to replicate the results obtained in the paper "Estimation of time-varying state correlations in state space models" by Caterina Schiavoni, Siem Jan Koopman, Franz Palm, Stephan Smeekes and Jan van den Brakel.


## R script *TV_corr.R*

This R script should be run in order to replicate the analyses, since it calls the *TV_corr.cpp* script. 

### Functions
The R script contains the following functions:
* ``Weights``: it creates the cubic splines weights.
* ``link``: it is a link function that bounds its argument between -1 and 1.
* ``KF_CC_known_corr``: it performs the Kalman filter estimation (and log-likelihood evaluation) with known values for the time-varying correlation.
* ``KF_CC_splines``: it performs the Kalman filter estimation (and log-likelihood evaluation) of the cubic splines model.
* ``KF_CC_const``: it performs the Kalman filter estimation (and log-likelihood evaluation) of the model with a time-constant correlation.
* ``IndInf_CC``: it performs the estimation of the static parameters of the nonlinear model by indirect inference.

### Inputs
Some inputs are common to many functions. Below you can find the list of all inputs:
* ``len``: the sample size.
* ``knots``:

## C++ script *TV_corr.cpp*

This C++ script contains a first set of functions that are needed in order to perform some mathematical operations. The functions that are instead related to the estimation of the time-varying correlation are:
* ``stratifiedResampling_rcpp``: it performs startified resampling.
* ``KF_CC_splines_rcpp``: it performs the Kalman filter estimation (and log-likelihood evaluation) of the cubic splines model.
* ``ucminf_rcpp_splines``: it maximises the the log-likelihood function evaluated in the ``KF_CC_splines_rcpp`` function.
* ``KF_t_rcpp``: it computes the prediction step of the Kalman filter.
* ``boot_filter_CC_rcpp``: it performs the Rao-blackwellised bootstrap filter estimation of the state vector of the nonlinear model.


In case of queries you can contact the corresponding author at c[dot]schiavoni[at]maastrichtuniversity[dot]nl.
