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

  real<lower=0> sigma;
  //real<lower=0> sigma2;
}

transformed parameters {

  vector[T] lp;
  lp = rep_vector(log_unif, T);
  for (s in 1:T) {
    for (t in 1:T) {
      lp[s] = lp[s] + normal_lpdf(y[t] | t < s ? x[t]*beta1+alpha1: x[t]*beta2+alpha2, sigma  );
    }
  }
  
  
}
model {
  
  alpha1 ~ normal(0,10);
  alpha2 ~normal(0,10);

  
  beta1 ~ normal(0,5);
  beta2 ~ normal(0,5);
  
  sigma ~ normal(0,5);
  //  sigma2 ~ normal(0,5);


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
          y_hat[t]=normal_rng(x[t]*beta1+alpha1,sigma);
          log_lik[t]=normal_lpdf(y[t]|x[t]*beta1+alpha1,sigma);
        }
      else
        {
          log_lik[t]=normal_lpdf(y[t]|x[t]*beta2+alpha2,sigma);
          y_hat[t]=normal_rng(x[t]*beta2+alpha2,sigma);
        }
    }

  
}
