data{
  #include "../input/y.stan"
}


transformed data {
  real log_unif;
  log_unif = -log(T);
}
parameters {
  real mu;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
}

transformed parameters {
  vector[T] lp;
  lp = rep_vector(log_unif, T);
  for (s in 1:T) {
    for (t in 1:T) {
      lp[s] = lp[s] + normal_lpdf(y[t] | mu , t<s ? sigma1:sigma2 );
    }
  }
}
model {

  mu ~ normal(0,0.1);
  sigma1 ~ normal(0,0.1);
  sigma2 ~ normal(0,0.1);
  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  s = categorical_logit_rng(lp);
  array[T] real log_lik;
  for( t in 1:T)
    {
      if( t < s)
        {
          log_lik[t]=normal_lpdf(y[t] | mu,sigma1);
        }
      else
        {
          log_lik[t]=normal_lpdf(y[t] | mu,sigma2);
        }
    }
}
