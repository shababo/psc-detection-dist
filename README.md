# psc-detection-dist
Distribution of code implementing the Bayesian PSC detection algorithm described in Merel et al, J. Neuro. Meths., 2016. 

## Implemenation Notes
-	The user should always input values of time in the units of seconds. However, the code converts everything to the unit of samples before running inference. This is the most convenient arrangement because then we only need to change the value of params.dt when the data has different sampling rates.
