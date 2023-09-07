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

  {
    vector[T + 1] lp_e;
    vector[T + 1] lp_l;
    lp_e[1] = 0;
    lp_l[1] = 0;

    for (t in 1:T) {
      lp_e[t + 1] = lp_e[t] + normal_lpdf(y[t] | x[t]*beta1+alpha1,sigma1);
      lp_l[t + 1] = lp_l[t] + normal_lpdf(y[t] | x[t]*beta2+alpha2,sigma2);
    }

    lp = rep_vector(log_unif + lp_l[T + 1], T)
      + head(lp_e, T) - head(lp_l, T);
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
