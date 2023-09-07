functions {
  #include "exp_lib.stan"
}

data{
##include "../input/y.stan"
  int T;
  int N;
  vector[T*N] y;
  array[N] int GT;
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
  sigma ~normal(0,1);
  tau ~ cauchy(0,Ps_tau);
  smooth~normal(0,Ps_smooth);


  for(n in 1:N)
    {
      int T=GT[n];
      vector[T] yy=segment(y,pos,T);
      array[T] real x;
      for(i in 1:T)
        {
          x[i]=i-1;
        }
      matrix[T,T] gp_obj = gp(x,tau,smooth,sigma);
      yy ~ multi_normal_cholesky(mu,gp_obj);
    }
  
}

generated quantities{
  array[N] vector[T] y_hat;

  for(n in 1:N)
    {
      int T=GT[n];
      vector[T] yy=segment(y,pos,T);
      array[T] real x;
      for(i in 1:T)
        {
          x[i]=i-1;
        }

      y_hat[i] = gp_pred_rng(x, yy, x, tau, smooth, sigma, delta);
    }

  
  for(i in 1:N)
    
#  matrix[N,N] gp_obj = gp(x,tau,smooth,sigma);#
#  y_hats=multi_normal_cholesky_rng(mu,gp_obj);
}
