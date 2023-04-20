data{
  int T;
  vector[T] y;
}

parameters{
  real<lower=0> sigma;
  real<lower=0> v;
}

model{

  v ~ normal(0,10);
  sigma ~normal(0,5);
  y ~ student_t(v,0,sigma);
  //x ~ v(mu,sigma);
}

generated quantities{
  array[T] real log_lik;
  array[T] real y_hat;
  for(t in 1:T)
    {
      log_lik[t]=student_t_lpdf(y[t]|v,0,sigma);
      y_hat[t]=student_t_rng(v,0,sigma);
    }
}
