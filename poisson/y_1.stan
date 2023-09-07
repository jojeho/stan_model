data{
  int T;
  vector[T] y;
}

parameters{
  real mu;
}

model{

  mu ~ normal(0,0.2);
  sigma ~normal(0,0.3);
  y ~ poisson(mu,sigma);
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

