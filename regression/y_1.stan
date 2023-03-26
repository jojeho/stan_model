data{
#include "../input/y.stan"
  vector[T] x;  
}

transformed data{
}

parameters{
  real beta;
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
      y[t]~ normal(x[t]*beta+alpha,sigma);
    }
}


/* generated quantities{ */
/*   vector[T] oblik; */
/*   vector[T] y_hat; */
/*   for(t in 1:T) */
/*     { */
/*       oblik[t]= normal_lpdf(y[t]|x[t]*beta+alpha,sigma); */
/*       y_hat[t]= normal_rng(x[t]*beta+alpha,sigma); */
/*     } */
/* } */
