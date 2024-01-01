data{
  #include "../input/y.stan"
}


transformed data {
  real log_unif;
  log_unif = -log(T);
}
parameters {
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  real<lower=0> v;
}

transformed parameters {

  vector[T] lp;
  {
    vector[T + 1] lp_e;
    vector[T + 1] lp_l;
    lp_e[1] = 0;
    lp_l[1] = 0;
    
    for (t in 1:T) {
      lp_e[t + 1] = lp_e[t] + student_t_lpdf(y[t] |v, 0,sigma1);
      lp_l[t + 1] = lp_l[t] + student_t_lpdf(y[t] |v, 0,sigma2);
    }
    lp = rep_vector(log_unif + lp_l[T + 1], T)
      + head(lp_e, T) - head(lp_l, T);
  }   
}
model {
  v  ~ normal(4,5);
  sigma1 ~ normal(0,1);
  sigma2 ~ normal(0,1);
  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  array[T] real log_lik;  
  s = categorical_logit_rng(lp);
}

