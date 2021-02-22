# Branson-2019-Regression-Discontinuity-Replication-Files

Included are the files needed to replicate the analysis in Branson et al. (2019), "A nonparametric Bayesian methodology for regression discontinuity designs."

Specifically, you should be able to replicate the analysis for Minutes Played (results shown in Table 1); the file pickMeans.mp.rds includes the data for this analysis. This also acts as an example so that you can see how to implement our methodology for other regression discontinuity designs.

All of the code is in R (branson2019RDDCode.R) and uses stan to fit the model (twoGPsOnlyCovarianceParamsWithSlope.stan). The comments in the R file should inform you what the functions are and how they operate. The end of the file (on line 170 of the code) demonstrates how to use the functions to estimate average treatment effects in a sharp RDD.

Note that you will need the rstan package to fit our model.

Of course, if you end up using our code, please cite our paper whenever you present results.
