# MNP: R package for fitting multinomial probit models using Markov chain Monte Carlo [![Build Status](https://travis-ci.org/kosukeimai/MNP.svg?branch=master)](https://travis-ci.org/kosukeimai/MNP) [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MNP)](https://cran.r-project.org/package=MNP) ![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/MNP)

This package fits the Bayesian multinomial probit model via Markov chain Monte Carlo.  The multinomial probit model is often used to analyze the discrete choices made by individuals recorded in survey data.  Examples where the multinomial probit model may be useful include the analysis of product choice by consumers in market research and the analysis of candidate or party choice by voters in electoral studies. The MNP package can also fit the model with different choice sets for each individual, and complete or partial individual choice orderings of the available alternatives from the choice set. The estimation is based on the efficient marginal data augmentation algorithm that is developed by Imai and van Dyk (2005). "[A Bayesian Analysis of the Multinomial Probit Model Using the Data Augmentation](https://doi.org/10.1016/j.jeconom.2004.02.002)," Journal of Econometrics, Vol. 124, No. 2 (February), pp. 311-334.  Detailed examples are given in Imai and van Dyk (2005). "[MNP: R Package for Fitting the Multinomial Probit Model](https://doi.org/10.18637/jss.v014.i03),"  Journal of Statistical Software, Vol. 14, No. 3 (May), pp. 1-32. 
