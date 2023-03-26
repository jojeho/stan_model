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
  int N;
  vector[N] y;
}

parameters{
  real mu;
  real<lower=0> sigma;
}

model{

  mu ~ normal(0,4);
  sigma ~normal(0,4);
  y ~ normal(mu,sigma);
  //x ~ v(mu,sigma);
}

generated quantities{
  array[N] real log_lik;
  array[N] real y_hat;
  for(t in 1:N)
    {
      log_lik[t]=normal_lpdf(y[t]|mu,sigma);
      y_hat[t]=normal_rng(mu,sigma);
    }
}
