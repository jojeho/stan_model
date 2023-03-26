data{

#include "../../input/X.stan"

}

parameters{
  
  vector[D] alpha;
  vector[D] beta;
  real<lower=0> sigma;
  real<lower=0> ps_beta;
  real<lower=0> ps_alpha;
  real alpha_z;
}

model{

  ps_beta ~ normal(0, 10);
  ps_alpha ~normal(0,10);
  alpha_z  ~ normal(0,1);

  for(d in 1:D)
    {
      beta[d] ~ normal(0,ps_beta);
      alpha[d] ~normal(0,ps_alpha);
    }

  #sigma ~normal(0,1);
  
  for(d in 1:D)
    {
      for(n in 2:N)
        X[n,d] ~ normal(X[n-1,d]*beta[d]+alpha[d] +alpha_z, sigma);
    }

}

  

