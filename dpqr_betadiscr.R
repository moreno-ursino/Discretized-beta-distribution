#######################################################################
#
# dbetadiscr(k,n,mu,sigma)
#
# k = outcome (between 0 and n)
# size = number of possibility
# mu = mean of the latent beta distribution
# sigma = precision parameter (alpha + beta) of the lantent
#	  beta distribution
######################################################################

dbetadiscr <<- function(x,size,mu,sigma,log=FALSE){

	if (sum(x<0)!=0 | sum(x>size)!=0 | sum(round(x))!= sum(x)) 
		cat ("ERROR: x must be integer between 0 and n \n")
	else
	if (size<=0 | round(size)!= size) 
		cat ("ERROR: size must be integer greater than 0 \n")
	else
	if (sum(mu<0)!=0 | sum(mu>1)!=0) 
		cat ("ERROR: mu must be between 0 and 1 \n")
	else
	if (sum(sigma<0)!=0 ) 
		cat ("ERROR: sigma must be positive \n")

	else
	pbeta((x+1)/(size+1),mu*sigma, (1-mu)*sigma) -
	pbeta((x)/(size+1),mu*sigma, (1-mu)*sigma)
}


#######################################################################
#
# pbetadiscr(x,size,mu,sigma)
#
# x = real number or vector of real
# size = number of possibility
# mu = mean of the latent beta distribution
# sigma = precision parameter (alpha + beta) of the lantent
#	  beta distribution
######################################################################

pbetadiscr <<- function(x,size,mu,sigma,lower.tail = TRUE, log.p = FALSE){

	if (size<=0 | round(size)!= size) 
		cat ("ERROR: size must be integer greater than 0 \n")
	else
	if (sum(mu<0)!=0 | sum(mu>1)!=0) 
		cat ("ERROR: mu must be between 0 and 1 \n")
	else
	if (sum(sigma<0)!=0 ) 
		cat ("ERROR: sigma must be positive \n")
	else {
	k = floor(x)
	pbeta((k+1)/(size+1),mu*sigma, (1-mu)*sigma,lower.tail = lower.tail, 
        log.p = log.p)		
 	}
}




#######################################################################
#
# qbetadiscr(u,n,mu,sigma)
#
# u = quantile (between 0 and 1)
# size = number of possibility
# mu = mean of the latent beta distribution
# sigma = precision parameter (alpha + beta) of the lantent
#	  beta distribution
######################################################################

qbetadiscr <<- function(p,size,mu,sigma,lower.tail = TRUE, log.p = FALSE){

	if (sum(p<0)!=0 | sum(p>1)!=0) 
		cat ("ERROR: p must be between 0 and 1 \n")
	else
	if (size<=0 | round(size)!= size) 
		cat ("ERROR: size must be integer greater than 0 \n")
	else
	if (sum(mu<0)!=0 | sum(mu>1)!=0) 
		cat ("ERROR: mu must be between 0 and 1 \n")
	else
	if (sum(sigma<0)!=0 ) 
		cat ("ERROR: sigma must be positive \n")
	else
	{
	ris = ceiling(qbeta(p, mu*sigma, (1-mu)*sigma,lower.tail = lower.tail, 
        log.p = log.p) * (size+1))-1
	ris[ris==-1] = rep(0,length(ris[ris==-1]))
	ris
	}
}



#######################################################################
#
# rbetadiscr(n,size,mu,sigma)
#
# n = number of observations. If length(n) > 1, the length is taken to be 
#	the number required.
# size = number of possibility
# mu = mean of the latent beta distribution
# sigma = precision parameter (alpha + beta) of the lantent
#	  beta distribution
######################################################################


rbetadiscr <<- function(n,size,mu,sigma){

	if (n<=0 | round(n)!= n) 
		cat ("ERROR: n must be integer greater than 0 \n")
	else
	if (size<=0 | round(size)!= size) 
		cat ("ERROR: size must be integer greater than 0 \n")
	else
	if (mu<0 | mu>1) 
		cat ("ERROR: mu must be between 0 and 1 \n")
	else
	if (sigma<0) 
		cat ("ERROR: sigma must be positive \n")
	else
	{
	a <- rbeta(n, mu*sigma, (1-mu)*sigma)
	sol <- floor(a*(size+1))
	sol[sol==size+1] = rep(size,length(sol[sol==size+1]))		
	sol
	}
}

