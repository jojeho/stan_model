data{
  #include "../../input/y.stan"
}

transformed data {
  real log_unif;
  vector[T] x;
  for(t in 1:T)
    {
      x[t]=t*0.1;
    }
    
  log_unif = -log(T);
}
parameters {
  real beta1;
  //real alpha1;
  
  real<lower=0> sigma;

  real<lower=0> kappa;
  real<lower=0> a;
  real<lower=0> b;  
}

transformed parameters {
  vector[T] lp;
  lp = rep_vector(log_unif, T);
  for (s in 1:T) {
    for (t in 1:T) {
      lp[s] = lp[s] + normal_lpdf(y[t] | t < s ? x[t]*beta1 : 0 , sigma);
    }
  }

}
model {

  //  alpha1 ~ normal(0,1);
  sigma ~ normal(0,3);

  kappa~gamma(a,b);
  a ~ gamma(2,1);
  b ~ gamma(2,1);
  beta1 ~ von_mises(0,kappa);
  
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
          log_lik[t]=normal_lpdf(y[t]|x[t]*beta1,sigma);
        }
      else
        {
          log_lik[t]=normal_lpdf(y[t]| 0,sigma);
        }
    }
}


