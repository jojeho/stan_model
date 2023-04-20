data{
  #include "../input/y.stan"
  real py;
}


transformed data {
  real log_unif;
  log_unif = -log(T);
}
parameters {
  real alpha1;
  real alpha2;
  real beta1;
  real beta2;
  real<lower=0> sigma;
}

transformed parameters {

  vector[T] lp;
  {
    vector[T + 1] lp_e;
    vector[T + 1] lp_l;
    lp_e[1] = 0;
    lp_l[1] = 0;

    lp_e[1 + 1] = lp_e[1] + normal_lpdf(y[1] | alpha1+beta1*py,sigma);
    lp_l[1 + 1] = lp_l[1] + normal_lpdf(y[1] | alpha2+beta2*py,sigma);
    
    for (t in 2:T) {
      lp_e[t + 1] = lp_e[t] + normal_lpdf(y[t] | alpha1+beta1*y[t-1],sigma);
      lp_l[t + 1] = lp_l[t] + normal_lpdf(y[t] | alpha2+beta2*y[t-1],sigma);
    }
    lp = rep_vector(log_unif + lp_l[T + 1], T)
      + head(lp_e, T) - head(lp_l, T);
  }   
}
model {
  
  alpha1 ~ normal(0,0.7);
  alpha2 ~ normal(0,0.7);
  sigma ~ normal(0,0.1);
  beta1 ~ normal(0,1);
  beta2 ~ normal(0,1);
  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  array[T] real log_lik;  
  s = categorical_logit_rng(lp);

  if(1 < s)
    log_lik[1]=normal_lpdf(y[1] | alpha1+beta1*py,sigma);
  else
    log_lik[1]=normal_lpdf(y[1] | alpha2+beta2*py,sigma);
  
  for( t in 2:T)
    {
      if( t < s)
        {
          log_lik[t]=normal_lpdf(y[t] | alpha1+beta1*y[t-1],sigma);
        }
      else
        {
          log_lik[t]=normal_lpdf(y[t] | alpha2+beta2*y[t-1],sigma);
        }
    }
  
}
