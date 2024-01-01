functions{
  void normal_lp(real mu)
  {
    mu ~ normal(0,10);
  }

  real v_lpdf(vector x , real mu,real sigma)
  {
    return normal_lpdf(x|mu,sigma);
  }
    
}

data{
  int T;
  vector[T] y;
}

parameters{
  real mu;
  real<lower=0> sigma;
}

model{

  mu ~ normal(0,0.03);
  sigma ~normal(0.3,0.05);
  y ~ normal(mu,sigma);
  //x ~ v(mu,sigma);
}

generated quantities{
  array[T] real log_lik;
  array[T] real y_hat;
  for(t in 1:T)
    {
      log_lik[t]=normal_lpdf(y[t]|mu,sigma);
      y_hat[t]=normal_rng(mu,sigma);
    }
}

