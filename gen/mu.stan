#include "input.lstan"

parameters{
  real mu;
  real<lower=0> sigma;
}

model{

  mu ~ normal(0,4);
  sigma ~normal(0,4);
  y ~ normal(mu,sigma);
  //x ~ v(mu,sigma);
}

generated quantities{
  array[NT] real log_lik;
  array[NT] real y_hat;
  for(t in 1:NT)
    {
      log_lik[t]=normal_lpdf(n_y[t]|mu,sigma);
      y_hat[t]=normal_rng(mu,sigma);
    }
}
