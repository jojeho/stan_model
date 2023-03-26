data{

#include "../../input/X.stan"
 array[N] real x;
}

transformed data{
  real lambda_mean =0.5;
  real period_scale=10;
  real  lambda_a = -lambda_mean * 2 / (lambda_mean - 1);
  real<lower=0> inv_period_scale = 1.0 / period_scale; 
}


parameters{
  //vector[D] beta;
  //vector<lower=0>[D] sigma;
  real sigma;
  // vector[D] alpha;

  real lambda;
  real rho;
  real theta_cycle;

  vector[D] kappa;
  vector[D] kappa_star;
}


transformed parameters{
  array[D] real    omega; // Cyclicality at time
  array[D] real    omega_star; // Cyclicality at time
  real rho_cos_lambda = rho*cos(lambda);
  real rho_sin_lambda = rho*sin(lambda);
  omega[1]=kappa[1];
  omega_star[1]=kappa_star[1];
  for(t in 2:D)
    {
      omega[t] = (rho_cos_lambda .* omega[t-1]) + (rho_sin_lambda .* omega_star[t-1]) + kappa[t];
      omega_star[t] = -(rho_sin_lambda .* omega[t-1]) + (rho_cos_lambda .* omega_star[t-1]) + kappa_star[t];
    }
}

model{

  rho ~ normal(0, 4);

  theta_cycle ~ cauchy(0, 2);
  
  (lambda / pi()) ~ beta(0.02, 2);
  kappa ~ normal(0, theta_cycle);
  kappa_star ~ normal(0, theta_cycle);
  
  
  sigma ~ cauchy(0,5);
  //alpha ~ normal(0,5);
  //beta  ~ normal(0,5);

  for(d in 1:D)
    {
      for(n in 1:N)
        {
          X[n,d]~ normal(x[n]*omega[d],sigma);
        }
    }
}

