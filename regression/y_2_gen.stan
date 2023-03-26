data{
#include "../input/y.stan"
  int D;
 array[T] vector[D] x;  
}

transformed data{
}

parameters{
  row_vector[D] beta;
  real<lower=0> sigma;
  //  vector[D] alpha;
  real alpha;
}

model{
  sigma ~ cauchy(0,5);
  alpha ~ normal(0,5);
  //  mu ~ normal(0,5);
  beta ~ normal(0,5);

  for(t in 1:T)
    {
      y[t]~ normal(beta*x[t]+alpha,sigma);
    }
}


generated quantities{
  vector[T] log_lik;
  vector[T] y_hat;
  for(t in 1:T)
    {
      log_lik[t]= normal_lpdf(y[t]|beta*x[t]+alpha,sigma);
      y_hat[t]= normal_rng(beta*x[t]+alpha,sigma);
    }
}
        
