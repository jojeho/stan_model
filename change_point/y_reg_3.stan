data{
  #include "../input/y.stan"
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

  real scale;
}

transformed parameters {

}
model {

  alpha1 ~ normal(0,0.1);
  alpha2 ~ normal(0,0.1);

  beta1 ~ normal(0,1);
  beta2 ~ normal(0,1);
  
  sigma1 ~ normal(0,1);
  sigma2 ~ normal(0,1);

  scale ~ normal(0,1);
  
  vector[T] lp;
  vector[T] new_y;
  for(t in 1:T)
    {
      new_y[t]=y[t]*scale;
    }
      
  lp = rep_vector(log_unif, T);
  for (s in 1:T) {
    for (t in 1:T) {
      lp[s] = lp[s] + normal_lpdf(new_y[t] | t < s ? x[t]*beta1+alpha1 : x[t]*beta2+alpha2 , t<s ? sigma1:sigma2 );
    }
  }
  
  target += log_sum_exp(lp);
}

/* generated quantities { */
/*   int<lower=1, upper=T> s; */
/*   s = categorical_logit_rng(lp); */
/* } */
