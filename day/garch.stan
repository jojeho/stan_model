data{
  int N;
  int R;
  matrix[N,R] X;

}

parameters{
  real mu ;
  //  real<lower=0.01> sigma;
  real alpha;
  vector[N] y;
  real<lower=0> sigma1;

  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=0> ss;
 }

transformed parameters {
  real<lower=0> sigma[N];
  
  sigma[1] = sigma1;
  for (t in 2:N)
    sigma[t] = sqrt(alpha0
                    + alpha1 * pow(y[t-1] - mu, 2)
                    + beta1 * pow(sigma[t-1], 2));  

}


model{

  alpha0~normal(0,1);
  mu ~normal(0,1);
  sigma ~ normal(0,2);
  sigma1 ~normal(0,1);
  alpha ~ normal(0,10);
  y ~ normal(0,1);
  ss ~ normal(0,1);
  for(n in 1:N)
    {
      X[n,] ~ normal(y[n],ss);
    }
}

