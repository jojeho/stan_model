data{
  #include "../../input/y.stan"
}

transformed data {
  real log_unif;
  vector[T] x;
  for(t in 1:T)
    {
      x[t]=t*0.1;
    }
  log_unif = -log(T);
  
}
parameters {
  real beta2;
  real beta1;
  //  real alpha2;
  real<lower=0> sigma1;
  real<lower=0> sigma2;

  real<lower=0> kappa;
  //  real<lower=0> a;
  //real<lower=0> b;
  
}


transformed parameters {

  vector[T] lp;
  lp = rep_vector(log_unif, T);
  for (s in 1:T) {
    for (t in 1:T) {
      lp[s] = lp[s] + normal_lpdf(y[t] | t < s ? x[t]*beta1  : x[t]*beta2  ,t<s ? sigma1: sigma2);
    }
  }
  
}

model {

  //alpha2 ~ normal(0,1);
  sigma1 ~ normal(0,1);
  sigma2 ~ normal(0,1);

  beta1 ~ von_mises(-pi()/4,kappa);
  beta2 ~ von_mises(pi()/4,kappa);
  kappa~normal(8,1);
  //kappa2~normal(2,2);
  //a ~ gamma(2,1);
  //b ~ gamma(2,1);
  
  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  array[T]  real log_lik;
  s = categorical_logit_rng(lp);
  for(t in 1:T)
    {
      if( t<s)
        {
          log_lik[t]=normal_lpdf(y[t]| x[t]*beta1,sigma1);
        }
      else
        {
          log_lik[t]=normal_lpdf(y[t]|x[t]*beta2,sigma2);
        }
    }
}
