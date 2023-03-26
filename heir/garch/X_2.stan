data{

#include "../../input/X.stan"

}


parameters{
  real mu[D];
  vector<lower=0>[D] alpha0;
  vector<lower=0,upper=1>[D] alpha1;
  vector<lower=0,upper=(1-alpha1)>[D] beta1;
  vector<lower=0>[D] sigma1;

  real mu_ps;
  real<lower=0> alpha0_ps;
  real<lower=0,upper=1> alpha1_ps;
  real<lower=0,upper=1> beta1_ps;
  real<lower=0> sigma1_ps;

  real mu_c;
  
}

transformed parameters{

}

model{
  
  alpha0_ps ~cauchy(0,5);
  alpha1_ps ~cauchy(0,5);
  beta1_ps ~cauchy(0,5);
  sigma1_ps ~cauchy(0,5);
  mu_ps~normal(0,1);
  mu_c ~normal(0,1);

  for(d in 1:D)
    {
      alpha0[d] ~ normal(0,alpha0_ps);
      alpha1[d] ~ normal(0,alpha1_ps);
      beta1[d] ~ normal(0,beta1_ps);
      sigma1[d] ~ normal(0,sigma1_ps);
      mu[d] ~ normal(0,mu_ps);
    }


  
  for(d in 1:D)
    {
      vector[N] sigma;
      sigma[1]=sigma1[d];;
      for(n in 2:N)
        sigma[n] = sqrt(alpha0[d]
                        + alpha1[d] * pow(X[n-1,d] - mu[d]+mu_c, 2)
                        + beta1[d] * pow(sigma[n-1], 2));
      X[,d] ~ normal(mu[d],sigma);
    }

}

  

