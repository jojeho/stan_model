

data{
#include "../../input/X.stan"

}


parameters{
  real<lower=0> alpha0;
  real alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;

  //real mu_ps;
  /* real<lower=0> alpha0_ps; */
  /* real<lower=0,upper=1> alpha1_ps; */
  /* real<lower=0,upper=1> beta1_ps; */
  //real<lower=0> sigma1_ps;

  vector[D] mu;
  /* real mu_c; */
  /* vector[D] mu_tilde; */
  /* real mu_tau; */
  vector<lower=0>[D] sigma1;
  /* real sigma_c; */
  /* vector[D] sigma_tilde; */
  /* real sigma_tau;   */
}

transformed parameters{


    
  /* for(d in 1:D) */
  /*   { */
  /*     mu[d]=mu_c + mu_tilde[d]*mu_tau; */
  /*     sigma1[d]=sigma_c + sigma_tilde[d]*sigma_tau; */
  /*   } */
}

model{
  
  /* alpha0_ps ~cauchy(0,5); */
  /* alpha1_ps ~cauchy(0,5); */
  /* beta1_ps ~cauchy(0,5); */
  //  sigma1_ps ~cauchy(0,5); 
  //mu_ps~normal(0,5);
  alpha1 ~ normal(0,5);
  beta1 ~ normal(0,5);
  alpha0 ~ normal(0,5);
  /* sigma_c ~ normal(0,2); */
  /* sigma_tau ~ cauchy(0,10); */
  /* mu_c ~ normal(0,1); */
  /* mu_tau ~ cauchy(0,5); */
  for(d in 1:D)
    {
      mu[d] ~ normal(0,1);
      sigma1[d]~normal(0,1);
    }

  
  for(d in 1:D)
    {
      vector[N] sigma;
      sigma[1]=sigma1[d];
      for(n in 2:N)
        sigma[n] = sqrt(alpha0
                        + alpha1 * pow(X[n-1,d] - mu[d], 2)
                        + beta1 * pow(sigma[n-1], 2));
      X[,d] ~ normal(mu[d],sigma);
    }

}

  

