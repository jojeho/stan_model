#include "input.lstan"

parameters {

  real mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=0> sigma1;
  
}

transformed parameters {
  
  real<lower=0> sigma[T];
  sigma[1] = sigma1;
  for (t in 2:T)
    sigma[t] = sqrt(alpha0
                    + alpha1 * pow(y[t-1] - mu, 2)
                    + beta1 * pow(sigma[t-1], 2));
}

model {
  alpha0~normal(0,1);
  
  mu ~normal(0,4);
  y~ normal(mu, sigma);
  sigma1 ~normal(0,3);
}

generated quantities{
  array[NT] real y_hat;
  array[NT] real log_lik;
  real<lower=0> n_sigma[NT];
  n_sigma[1] = sigma1;
  for (t in 2:NT)
    n_sigma[t] = sqrt(alpha0
                    + alpha1 * pow(n_y[t-1] - mu, 2)
                    + beta1 * pow(n_sigma[t-1], 2));

  
  for(t in 1:NT)
    {
      y_hat[t]=normal_rng(mu,n_sigma[t]);
      log_lik[t]=normal_lpdf(n_y[t]|mu,n_sigma[t]);
    }
  //  y_hat=sigma;
}

