data{
  #include "../../input/y.stan"
}

transformed data {
  real log_unif;
  vector[T] x;
  for(t in 1:T)
    {
      x[t]=t*0.01;
    }
  log_unif = -log(T);
  
}
parameters {
  real beta2;
  //  real alpha2;
  real<lower=0> sigma;
  
}

transformed parameters {

  vector[T] lp;
  lp = rep_vector(log_unif, T);
  for (s in 1:T) {
    for (t in 1:T) {
      lp[s] = lp[s] + normal_lpdf(y[t] | t < s ? 0  : x[t]*beta2  , sigma );
    }
  }

  
}

model {

  //alpha2 ~ normal(0,1);
  beta2 ~ normal(0,1);
  sigma ~ normal(0,3);


  
  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  array[T]  real log_lik;
  s = categorical_logit_rng(lp);
  for(t in 1:T)
    {
      if( t<s)
        {
          log_lik[t]=normal_lpdf(y[t]| 0,sigma);
        }
      else
        {
          log_lik[t]=normal_lpdf(y[t]|x[t]*beta2,sigma);
        }
    }
}
