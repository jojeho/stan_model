data {
  int T;
  array[T] real y;

  int NT;
  array[NT] real ny;
  //real PM_alpha;
  //real PS_alpha;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(0,1);
  beta ~ normal(0,1);
  for (n in 2:T) {
    y[n] ~ normal(alpha + beta * y[n-1], sigma);
  }
}



generated quantities{
  array[NT-1] real log_lik;
  array[NT-1] real y_hat;
  for(n in 2:NT)
    {
      log_lik[n-1]=normal_lpdf(ny[n]|alpha + beta*ny[n-1],sigma);
      y_hat[n-1]=normal_rng(alpha + beta*ny[n-1],sigma);
    }  
}
