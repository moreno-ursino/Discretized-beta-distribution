#################################################################
#
#		Script written by Ursino 5/11/2014
#
#
# Dicrete-beta with 1 random effect: estimation using Stan
# 
##################################################################


library(rstan)
sm3 <- stan_model(file="dbeta_rando.stan", model_name='p2mod',verbose=FALSE)

source("dpqr_betadiscr.R")
source("discrete_beta_R.R")


# example

data_ex = data.frame(
  y = c(4,5,6,7,5,6,7,4,3,4,5,6,5,9,4,8,3,7,6,9,3,6,2,6),
  time = factor(rep(1:2, 12)),
  group = factor(rep(1:2, each=12)),
  id = rep(1:12, each=2)
)

y = cbind(data_ex$y,10-data_ex$y)  # assuming 10 as maximum value for y 
formulamu = ~  group + time 
formulaphi = ~ time
id = as.integer(data_ex$id)


size <- sum(y[1,]) # size of the betabinomial

#model matrices 
mat_mutot <- model.matrix(formulamu, data_ex)
mat_phitot <- model.matrix(formulaphi, data_ex)
K1=dim(mat_mutot)[2] 
K2=dim(mat_phitot)[2]       
n = length(unique(data_ex$id))

# using empirical Bayes to find prior means
starting <- bdmod(y,formulamu,formulaphi,data_ex)

# data for Stan

# N // number of patient
# K1 // number of covariate for mu
# K2 // number of covariate for phi
# t // number of times
# Size // maximum outcome vale on a range from 0 to Size
# matmu // model matrix for mu
# matphi // model matrix for phi
# patient[N*t] // id
# y[N*t] // outcomes
# smu[K1] // prior means for mu covariates
# sphi[K2] // prior means for \phi covariates


data_s <- list(N=n,K1=K1, K2=K2, t=2, Size=size, matmu=mat_mutot, matphi=mat_phitot, 
               patient=id, y=y[,1], smu=starting$MU[,1], sphi=starting$PHI[,1] )

reg <- sampling(sm3, data=data_s,iter=2000, chains=2)

# Suggestions:
# check convergence (iter=2000 and 2 chains usually are not enough), 
# change priors definition in the dbeta_rando.stan file, 
# change parametrizations, etc...




