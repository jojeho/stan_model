data{

#include "../../input/X.stan"

}

parameters{
  
  vector<lower=0>[D] sigma;
  real beta_c;
  real alpha_c;
  vector[D] beta_tilde;
  vector[D] alpha_tilde;
  real tau;
  real tau2;

  real alpha_z_tilde;
  real alpha_z_c;
  real alpha_z_tau;
}

transformed parameters{
  vector[D] beta;
  vector[D] alpha;
  real alpha_z;
  for(d in 1:D)
    {
      beta[d]=beta_c+beta_tilde[d]*tau;
      alpha[d]=alpha_c + alpha_tilde[d]*tau2;
    }

  alpha_z=alpha_z_c + alpha_z_tilde*alpha_z_tau;
}

model{
  tau ~ cauchy(0,10);
  tau2 ~ cauchy(0,10);
  #alpha_ps ~normal(0,5);
  alpha_z  ~ normal(0,1);
  beta_c ~ normal(1,10);
  alpha_c ~ normal(0,1);
  sigma ~normal(0,10);
  beta_tilde~normal(0,1);
  alpha_tilde~normal(0,1);

  alpha_z_tilde~normal(0,1);
  alpha_z_tau~ cauchy(0,5);
  alpha_z_c ~normal(0,2);
  
  for(d in 1:D)
    {
      for(n in 2:N)
        X[n,d] ~ normal(X[n-1,d]*beta[d]+alpha[d] +alpha_z, sigma[d]);
    }

}

  

