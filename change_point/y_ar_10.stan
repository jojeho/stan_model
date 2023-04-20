data{
  #include "../input/y.stan"
  real py;
}

transformed data {
  real log_unif;
  log_unif = -log(T);
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
}

transformed parameters {
  real alpha1;
  real alpha2;
  alpha1=alpha;
  alpha2=-alpha1;
    
  vector[T] lp;
  {
    vector[T + 1] lp_e;
    vector[T + 1] lp_l;
    lp_e[1] = 0;
    lp_l[1] = 0;

    lp_e[1 + 1] = lp_e[1] + normal_lpdf(y[1] | alpha1+beta*py,sigma1);
    lp_l[1 + 1] = lp_l[1] + normal_lpdf(y[1] | alpha2+beta*py,sigma2);
    
    for (t in 2:T) {
      lp_e[t + 1] = lp_e[t] + normal_lpdf(y[t] | alpha1+beta*y[t-1],sigma1);
      lp_l[t + 1] = lp_l[t] + normal_lpdf(y[t] | alpha2+beta*y[t-1],sigma2);
    }
    lp = rep_vector(log_unif + lp_l[T + 1], T)
      + head(lp_e, T) - head(lp_l, T);
  }   
}
model {
  alpha ~ normal(0,0.2);
  sigma1 ~ normal(0,1);
  sigma2 ~ normal(0,1);
  beta ~ normal(0,1);
  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  s = categorical_logit_rng(lp);
}
