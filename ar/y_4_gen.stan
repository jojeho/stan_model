data{

#include "../input/y.stan"
  array[T] real x;
}


parameters{
  
  real<lower=0> sigma;
  //real beta;
  real beta1;
  //  real gamma;
  //real alpha_tilde;
  //real tau;
  real alpha;
  real beta;
}

transformed parameters{

}


model{

  //alpha_tilde ~ normal(0,1);
 // tau ~ cauchy(0,5);
  alpha ~ normal(0,1);
  #beta~normal(0,2);
  beta1~normal(0,2);
  beta~normal(0,2);
  //beta1~gamma(2,2);
  sigma ~ cauchy(0,5);
  //  gamma~normal(0,1);
  
  array[T] real y_hat;
  //  real beta=exp(beta_log);
  for(t in 1:T)
    {
      y_hat[t] = beta1*exp(-x[t]*beta)+alpha;
    }
  
  for(t in 1:T)
    {
      y[t] ~ normal(y_hat[t],sigma);
    }
}

generated quantities{
  vector[T] oblik;
  vector[T] y_hat;
  //  real beta=exp(beta_log);
  for(t in 1:T)
    {
      oblik[t] =normal_lpdf(y[t]| beta1*exp(-x[t]*beta)+alpha,sigma);
      y_hat[t] =normal_rng( beta1*exp(-x[t]*beta)+alpha,sigma);
    }
}
