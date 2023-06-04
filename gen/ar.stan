#include "input.lstan"

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
  real y_init;
}
model {
  sigma ~normal(0,3);
  alpha ~ normal(0,2);
  beta ~ normal(0,1);
  y_init ~ normal(0,2);
  y[1] ~ normal(alpha + beta*y_init,sigma);
  for (n in 2:T) {
    y[n] ~ normal(alpha + beta * y[n-1], sigma);
  }
}


generated quantities{
  array[NT] real log_lik;
  array[NT] real y_hat;

  log_lik[1]=normal_lpdf(n_y[1]|alpha + beta*y_init,sigma);
  y_hat[1]=normal_rng(alpha + beta*y_init,sigma);
  
  for(n in 2:NT)
    {
      log_lik[n]=normal_lpdf(n_y[n]|alpha + beta*n_y[n-1],sigma);
      y_hat[n]=normal_rng(alpha + beta*n_y[n-1],sigma);
    }  
}
