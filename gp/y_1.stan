functions {
  #include "exp_lib.stan"
}

data{
##include "../input/y.stan"
  int T;
  vector[T] y;
  real Ps_tau;
  real Ps_smooth;
}

transformed data{
  vector[T] mu;
  array[T] real x;
  for(i in 1:T)
    {
      x[i]=i-1;
    }
  
  real delta=1e-10;
  mu =rep_vector(0,T);
}

parameters{
  real<lower=0> tau;
  real<lower=0> smooth;
  real<lower=0> sigma;
}


model{
  sigma ~normal(0,0.1);
  tau ~ cauchy(0,Ps_tau);
  smooth~normal(0,Ps_smooth);

  matrix[T,T] gp_obj = gp(x,tau,smooth,sigma);
  y ~ multi_normal_cholesky(mu,gp_obj);
  
}

generated quantities{
  vector[T] y_hat;
  y_hat = gp_pred_rng(x, y, x, tau, smooth, sigma, delta);
#  matrix[N,N] gp_obj = gp(x,tau,smooth,sigma);#
#  y_hats=multi_normal_cholesky_rng(mu,gp_obj);
}
