data{
#include "../input/y.stan"
  vector[T] x;  
}

transformed data{
}

parameters{
  real beta;
  real sigma;
}

model{
  sigma ~ cauchy(0,5);
  beta ~ normal(0,5);

  for(t in 1:T)
    {
      y[t]~ normal(x[t]*beta,sigma);
    }
}


generated quantities{
  vector[T] oblik;
  vector[T] y_hat;
  for(t in 1:T)
    {
      oblik[t]= normal_lpdf(y[t]|x[t]*beta,sigma);
      y_hat[t]= normal_rng(x[t]*beta,sigma);
    }
}
