#################################################################
#
#		Script written by Ursino 5/11/2014
#
################################################################

# INPUT:
#	y		              ordinal variable (0,1,...m). Use cbind(variable, m - variable)
#	formulamu	        formula for \mu  "~ gender + age + ...". Do not add anything before ~
#	formulaphi	      formula for \phi  "~ gender + age + ...". Do not add anything before ~
#	data		          dataframe
#
#
# OUTPUT:
#	$MU		            estimates and p-value for \mu covariates
#	$PHI		          estimates and p-value for \phi covariates
#	AIC
#	AICc

# Example of use:
#		formulamu = ~ Gender + Age
#		formulaphi = ~ Gender
#		bdmod(cbind(y,10-y), formulamu, formulaphi, data)


##################################################################
library(numDeriv)


bdmod <- function(y ,formulamu, formulaphi, data)
{ 
  #design matrix for \mu
  mat_mu <- model.matrix(formulamu, data)
  mat_mup <<- mat_mu
  
  #design matrix for \phi
  mat_phi <- model.matrix(formulaphi, data)
  mat_phip <<- mat_phi
  
  n <- dim(data)[1] #number of individuals
  
  size <- sum(y[1,]) #n for the binomial likelihood
  
  
  # loglikelihood function
  
  logverosim <- function(par,x1,size,mat_mu,mat_phi){
    x=x1
    beta_mu <- par[1:dim(mat_mu)[2]]
    beta_phi <- par[(dim(mat_mu)[2]+1) : (dim(mat_mu)[2]+dim(mat_phi)[2])]
    
    # vectors are considered as column ones.
    
    mu <- 1/(1+exp(- (mat_mu %*% beta_mu))) #\mu vector: one for each individual
    phi <- exp(mat_phi %*% beta_phi) #\phi vector: one for each individual
    
    # numerical checking
    mu[mu==1] <- rep(1-.Machine$double.eps, length(mu[mu==1])) 
    mu[mu==0] <- rep(2^(-400), length(mu[mu==0])) 
    phi[phi==0] <- rep(2^(-400), length(phi[phi==0]))
    phi[phi==Inf] <- rep(2^(500), length(phi[phi==Inf]))
    
    alfa1 <- mu*phi
    alfa2 <- (1-mu)*phi
    
    # numerical checking
    alfa1[alfa1==0] <- rep(2^(-400), length(alfa1[alfa1==0]))
    alfa1[alfa1==Inf] <- rep(2^(500), length(alfa1[alfa1==Inf]))
    alfa2[alfa2==0] <- rep(2^(-400), length(alfa2[alfa2==0]))
    alfa2[alfa2==Inf] <- rep(2^(500), length(alfa2[alfa2==Inf]))
    
    
    a <- pbeta(x/(size +1),alfa1,alfa2) 
    b <- pbeta((x+1)/(size +1),alfa1,alfa2)
    l <- sum(log(b-a))
    -l
  }
  
  # estimation step
  stima = optim(par=rep(0, dim(mat_mu)[2] + dim(mat_phi)[2]) , fn=logverosim, x1=y[,1], size=size,
                mat_mu=mat_mu, mat_phi=mat_phi,  method = "BFGS", hessian=TRUE )
  
  # extracting results
  per_mu = matrix(stima$par[1:dim(mat_mu)[2]], ncol=1)
  per_mu = data.frame(per_mu)
  dimnames(per_mu)[[1]] <- dimnames(mat_mu)[[2]]
  dimnames(per_mu)[[2]] <- c("Estimate")
  
  per_phi = matrix(stima$par[(dim(mat_mu)[2]+1) : (dim(mat_mu)[2]+dim(mat_phi)[2])], ncol=1)
  per_phi = data.frame(per_phi)
  dimnames(per_phi)[[1]] <- dimnames(mat_phi)[[2]]
  dimnames(per_phi)[[2]] <- c("Estimate")
  
  
  # p-values
  simm <- stima$hessian
  
  simm <- hessian(func=logverosim, x=stima$par, x1=y[,1], size=size,
          mat_mu=mat_mu, mat_phi=mat_phi)
  
  simm <- solve(simm)
  sdval <- diag(simm)
  sdval <- sqrt(sdval)
  
  pvalmu <- 2*(1-pnorm(abs(per_mu[,1]/sdval[1:dim(per_mu)[1]]))) 
  pvalphi <- 2*(1-pnorm(abs(per_phi[,1]/sdval[-(1:dim(per_mu)[1])]))) 
  
  codemu <- rep(" ", dim(per_mu)[1])
  codemu[which(pvalmu<=0.001)] <- rep("***",length(pvalmu[pvalmu<=0.001]))
  codemu[which(pvalmu>0.001 & pvalmu<=0.01)] <- rep("**",length(which(pvalmu>0.001 & pvalmu<=0.01)))
  codemu[which(pvalmu>0.01 & pvalmu<=0.05)] <- rep("*",length(which(pvalmu>0.01 & pvalmu<=0.05)))
  codemu[which(pvalmu>0.05 & pvalmu<=0.1)] <- rep(".",length(which(pvalmu>0.05 & pvalmu<=0.1)))
  
  codephi <- rep(" ", dim(per_phi)[1])
  codephi[which(pvalphi<=0.001)] <- rep("***",length(pvalphi[pvalphi<=0.001]))
  codephi[which(pvalphi>0.001 & pvalphi<=0.01)] <- rep("**",length(which(pvalphi>0.001 & pvalphi<=0.01)))
  codephi[which(pvalphi>0.01 & pvalphi<=0.05)] <- rep("*",length(which(pvalphi>0.01 & pvalphi<=0.05)))
  codephi[which(pvalphi>0.05 & pvalphi<=0.1)] <- rep(".",length(which(pvalphi>0.05 & pvalphi<=0.1)))
  
  # final tables
  per_mu <- cbind(per_mu, Sd_error= sdval[1:dim(per_mu)[1]], 
                  Wald_Test=per_mu[,1]/sdval[1:dim(per_mu)[1]], Pval= pvalmu, codemu)
  per_phi <- cbind(per_phi, Sd_error=sdval[-(1:dim(per_mu)[1])],
                   Wald_Test=per_phi[,1]/sdval[-(1:dim(per_mu)[1])], Pval=pvalphi, codephi)
  
  dimnames(per_mu)[[2]][5] = " "
  dimnames(per_phi)[[2]][5] = " "
  
  k= length(stima$par)
  N <- dim(data)[1]
  AIC = 2*( k + stima$value) 
  AICc = AIC + 2*k*(k+1)/(N -k -1) 
  
  list( MU = per_mu, PHI = per_phi, AICc = AICc, AIC= AIC)
}
