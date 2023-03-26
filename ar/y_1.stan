
data {
  #include "../input/y.stan"
  real PM_alpha;
  real PS_alpha;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(PM_alpha,PS_alpha);
  beta ~ normal(0,1);
  for (n in 2:T) {
    y[n] ~ normal(alpha + beta * y[n-1], sigma);
  }
}



generated quantities{
  array[T] real log_lik;
  array[T] real y_hat;
  log_lik[1]=normal_lpdf(y[1]|y[1],sigma);
  y_hat[1]=normal_rng(y[1],sigma);
  for(n in 2:T)
    {
      log_lik[n]=normal_lpdf(y[n]|alpha + beta*y[n-1],sigma);
      y_hat[n]=normal_rng(alpha + beta*y[n-1],sigma);

    }
  
}
