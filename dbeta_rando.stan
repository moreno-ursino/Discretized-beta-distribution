data {
  int<lower=0> N; // number of patient
  int<lower=0> K1; // covariate mu
  int<lower=0> K2; // covariate phi
  int<lower=0> t; // number of times
  int<lower=0> Size; // size
  matrix[N*t,K1] matmu; // matmu
  matrix[N*t,K2] matphi; // matmu
  int patient[N*t];
  int y[N*t]; // outcomes
  real smu[K1]; // starting points for mu
  real sphi[K2]; // starting points for phi
}
parameters {
  vector[K1] betamu;
  vector[K2] betaphi;
  vector[N] bmu;
  real<lower=0.0001> sigma;
}
model {
  real Size1; // size
  vector[N*t] bmutot;
  vector[N*t] mu;
  vector[N*t] phi;
  vector[N*t] alfa1;
  vector[N*t] alfa2;
  vector[Size] p_pat1;
  vector[N*t] a;
  vector[N*t] b;
  
  Size1 = Size/(1.0);
  
  //random effects for each patient
  for (k in 1:(N*t)){
    bmutot[k] = bmu[patient[k]];
  }
  
  mu = (matmu * betamu) + bmutot;
  phi = exp(matphi * betaphi);   // vectorial cooding 
  
  for (k in 1:(N*t)){
    mu[k] = 1/(1+exp(- mu[k])); // mu and phi vectors
    mu[k] = if_else(mu[k]==0, 2^(-900), mu[k]);
    mu[k] = if_else(mu[k]==1, 1-machine_precision(), mu[k]); // numerical control 
    phi[k] = if_else(is_inf(phi[k])==1, 2^(900), phi[k]);
    phi[k] = if_else(phi[k]==0,  2^(-900), phi[k]);
  }
  
  
  alfa1 = mu .* phi;
  alfa2 = (1-mu).*phi;    
  
  // numerical control
  for (k in 1:(N*t)){
    alfa1[k] = if_else(is_inf(alfa1[k])==1, 2^(900), alfa1[k]);
    alfa1[k] = if_else(alfa1[k]==0,  2^(-900), alfa1[k]);
    alfa2[k] = if_else(is_inf(alfa2[k])==1, 2^(900), alfa2[k]);
    alfa2[k] = if_else(alfa2[k]==0,  2^(-900), alfa2[k]);
  }
  
  
  // likelihood
  for (i in 1:(N*t) ){
    
    a[i] = beta_cdf(y[i]/(Size1 +1),alfa1[i],alfa2[i]);
    b[i] = beta_cdf((y[i]+1)/(Size1 +1),alfa1[i],alfa2[i]);
    target += log(b[i]-a[i]);
    
  }
  
  // prior distributions
  betamu[1] ~ normal(smu[1], 100);
  for (l in 2:K1) betamu[l] ~ normal(smu[l], 10);
  betaphi[1] ~ normal(sphi[1], 100);
  for (l in 2:K2) betaphi[l] ~ normal(sphi[l], 10);
  bmu ~ normal(0, sigma);
  sigma ~ inv_gamma(0.01, 0.01);
}
