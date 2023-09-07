data {
  int<lower=0> T;
  real y[T];
  int N;
  int n;
  //  real<lower=0> sigma1;
}

parameters {

  real mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  //vector<lower=0>[N] sigma1;
  real<lower=0> sigma1;
  
}


generated quantities{
  array[T] real y_hat;
  array[T] real oblik;
  real<lower=0> sigma[T];
  //sigma[1] = sigma1[n];
  sigma[1]=sigma1;
  for (t in 2:T)
    sigma[t] = sqrt(alpha0
                    + alpha1 * pow(y[t-1] - mu, 2)
                    + beta1 * pow(sigma[t-1], 2));

  
  for(t in 1:T)
    {
      y_hat[t]=normal_rng(mu,sigma[t]);
      oblik[t]=normal_lpdf(y[t]|mu,sigma[t]);
    }
  //  y_hat=sigma;
}

