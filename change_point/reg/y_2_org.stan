data{
  #include "../../input/y.stan"
  vector[T] x;
}


transformed data {
  real log_unif;
  log_unif = -log(T);
}
parameters {
  real beta1;
  real beta2;

  real alpha1;
  real alpha2;

  real<lower=0> sigma1;
  real<lower=0> sigma2;
}

transformed parameters {

  vector[T] lp;
  lp = rep_vector(log_unif, T);
  for (s in 1:T) {
    for (t in 1:T) {
      lp[s] = lp[s] + normal_lpdf(y[t] | t < s ? x[t]*beta1+alpha1: x[t]*beta2+alpha2,
                                  t<s ? sigma1:sigma2);
    }
  }
  
  
}
model {
  
  alpha1 ~ normal(0,1);
  alpha2 ~normal(0,1);

  
  beta1 ~ cauchy(0,2.5);
  beta2 ~ cauchy(0,2.5);
  
  sigma1 ~ normal(0,1);
  sigma2 ~ normal(0,1);


  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  vector[T] y_hat;
  vector[T] log_lik;
  s = categorical_logit_rng(lp);
  for(t in 1:T)
    {
      if( t<s)
        {
          y_hat[t]=normal_rng(x[t]*beta1+alpha1,sigma1);
          log_lik[t]=normal_lpdf(y[t]|x[t]*beta1+alpha1,sigma1);
        }
      else
        {
          log_lik[t]=normal_lpdf(y[t]|x[t]*beta2+alpha2,sigma2);
          y_hat[t]=normal_rng(x[t]*beta2+alpha2,sigma2);
        }
    }

  
}
