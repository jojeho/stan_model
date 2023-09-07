data {
  int<lower=0> T;
  real y[T];
  //  real<lower=0> sigma1;
}

parameters {

  real mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=0> sigma1;
  
}

transformed parameters {
  
  real<lower=0> sigma[T];
  
  sigma[1] = sigma1;
  for (t in 2:T)
    sigma[t] = sqrt(alpha0
                    + alpha1 * pow(y[t-1] - mu, 2)
                    + beta1 * pow(sigma[t-1], 2));
}

model {

  //alpha0~normal(0,1);
  //alpha1 ~ normal(0,1);
  //beta ~ normal(0,1);
  //mu ~normal(0,1);
  y~ normal(mu, sigma);
  sigma1 ~normal(0,0.3);
}

/* generated quantities{ */
/*   array[T] real y_hat; */
/*   array[T] real oblik; */
/*   for(t in 1:T) */
/*     { */
/*       y_hat[t]=normal_rng(mu,sigma[t]); */
/*       oblik[t]=normal_lpdf(y[t]|mu,sigma[t]); */
/*     } */
/*   //  y_hat=sigma; */
/* } */

