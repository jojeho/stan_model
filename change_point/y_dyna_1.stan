data{
  #include "../input/y.stan"
}


transformed data {
  real log_unif;
  log_unif = -log(T);
}
parameters {
  real mu1;
  real mu2;
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
      lp_e[t + 1] = lp_e[t] + normal_lpdf(y[t] | mu1,sigma1);
      lp_l[t + 1] = lp_l[t] + normal_lpdf(y[t] | mu2,sigma2);
    }
    lp = rep_vector(log_unif + lp_l[T + 1], T)
      + head(lp_e, T) - head(lp_l, T);
  }   
}
model {

  mu1 ~ normal(0,1);
  mu2 ~ normal(0,1);
  sigma1 ~ normal(0,1);
  sigma2 ~ normal(0,1);
  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  array[T] real log_lik;  
  s = categorical_logit_rng(lp);
  for( t in 1:T)
    {
      if( t < s)
        {
          log_lik[t]=normal_lpdf(y[t] | mu1,sigma1);
        }
      else
        {
          log_lik[t]=normal_lpdf(y[t] | mu2,sigma2);
        }
    }
  
}
