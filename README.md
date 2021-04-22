In this folder, you can find the code to run the fixed effect model - the discretized beta distribution - presented in “A new parsimonious model for ordinal longitudinal data with application to subjective evaluations of a gastrointestinal disease”, Statistical methods in medical research 27.5 (2018): 1376-1393.
Moreover, a hint on how to write the Bayesian model in Stan is given.

The scripts can be run in all operating system. R must be installed with the following libraries: 
   numDeriv, rstan.  

Files description:
- dpqr_betadiscr.R: the probability mass function, cumulative distribution function, quantile and random function of the discretized beta distribution;
- discrete_beta_R.R: function for parameters' estimation;
- dbeta_rando.stan: stan model in case of one random effect for \mu;
- discrete_beta_random_eff.R: example of script using the fixed effect model to determine the means of the priors distributions for the mixed-effect model;

Attention: dbeta_rando.stan must be in the same folder where you are running the scripts. To modify the prior distributions and/or the model parametrization, you need to work on this file for the mixed-effect model.  

For any issue or question, please contact moreno.ursino@gmail.com.
