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
  lp = rep_vector(log_unif, T);
  for (s in 1:T) {
    for (t in 1:T) {
      lp[s] = lp[s] + normal_lpdf(y[t] | t < s ? mu1 : mu2 , t<s ? sigma1:sigma2 );
    }
  }
}
model {

  mu1 ~ normal(0,0.1);
  mu2 ~ normal(0,0.1);
  sigma1 ~ normal(0,0.1);
  sigma2 ~ normal(0,0.1);
  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  s = categorical_logit_rng(lp);
}
