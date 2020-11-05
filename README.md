# Estimation of time-varying state correlations in state space models

This repository contains the R and C++ codes that can be used to replicate the results obtained in the paper "Estimation of time-varying state correlations in state space models" by Caterina Schiavoni, Siem Jan Koopman, Franz Palm, Stephan Smeekes and Jan van den Brakel.


## R script *TV_corr.R*

This R script should be run in order to replicate the analyses, since it calls the *TV_corr.cpp* script. It contains the following functions:
* *Weights*: it creates the cubic splines weights.
* *link*: it is a link function that bounds its argument between -1 and 1.
* *KF_CC_known_corr*: it performs the Kalman filter estimation (and log-likelihood evaluation) with known values for the time-varying correlation.
* *KF_CC_splines*: it performs the Kalman filter estimation (and log-likelihood evaluation) of the cubic splines model.
* *KF_CC_const*: it performs the Kalman filter estimation (and log-likelihood evaluation) of the model with a time-constant correlation.
* *IndInf_CC*: it performs the estimation of the static parameters of the nonlinear model by indirect inference.

## C++ script *TV_corr.cpp*


In case of queries you can contact the corresponding author at c[dot]schiavoni[at]maastrichtuniversity[dot]nl.
