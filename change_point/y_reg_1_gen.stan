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
}

transformed parameters {

}


model {

  alpha1 ~ normal(0,1);
  alpha2 ~ normal(0,1);

  beta1 ~ normal(0,1);
  beta2 ~ normal(0,1);
  
  sigma1 ~ normal(0,1);
  sigma2 ~ normal(0,1);

  vector[T] lp;
  lp = rep_vector(log_unif, T);
  for (s in 1:T) {
    for (t in 1:T) {
      lp[s] = lp[s] + normal_lpdf(y[t] | t < s ? x[t]*beta1+alpha1 : x[t]*beta2+alpha2 , t<s ? sigma1:sigma2 );
    }
  }
  
  target += log_sum_exp(lp);
}


generated quantities { 
  int<lower=1, upper=T> s;
  vector[T] log_lik;
  vector[T] y_hat;
  
  oblik = rep_vector(log_unif, T);
  for (m in 1:T) {
    for (t in 1:T) {
      log_lik[m] = log_lik[m] + normal_lpdf(y[t] | t < m ? x[t]*beta1+alpha1 : x[t]*beta2+alpha2 , t<m ? sigma1:sigma2 );
    }
  }

  s = categorical_logit_rng(oblik);
  for(t in 1:T)
    {
      if( t < s)
        {
          log_lik[t]=normal_lpdf(y[t]|x[t]*beta1+alpha1 ,sigma1);
          y_hat[t]=normal_rng(x[t]*beta1+alpha1 ,sigma1);
        }
      else
        {
          log_lik[t]=normal_lpdf(y[t]| x[t]*beta2+alpha2 , sigma2);
          y_hat[t]=normal_rng( x[t]*beta2+alpha2 , sigma2);
        }
    }

    
  
} 
