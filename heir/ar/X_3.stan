data{

#include "../../input/X.stan"

}

parameters{
  
  vector[D] alpha;
  
  vector<lower=0>[D] sigma;
  real<lower=0> alpha_ps;
  real alpha_z;
  real beta_c;
  vector[D] beta_tilde;
  real tau;
  
}

transformed parameters{
  vector[D] beta;
  for(d in 1:D)
    {
      beta[d]=beta_c+beta_tilde[d]*tau;
    }
}

model{
  tau ~ cauchy(0,10);
  alpha_ps ~normal(0,5);
  alpha_z  ~ normal(0,1);
  beta_c ~ normal(0,1);
  for(d in 1:D)
    {
      alpha[d] ~normal(0,alpha_ps);
    }

  sigma ~normal(0,10);
  
  for(d in 1:D)
    {
      for(n in 2:N)
        X[n,d] ~ normal(X[n-1,d]*beta[d]+alpha[d] +alpha_z, sigma[d]);
    }

}

  

