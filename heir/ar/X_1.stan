data{

#include "../../input/X.stan"

}

parameters{
  
  //  vector[D] alpha;
  vector[D] beta;
  vector<lower=0>[D] sigma;

  /* real mu_c; */
  /* real mu_tilde; */
  /* real mu_tau; */
  real alpha;
}

transformed parameters{
   /* real mu; */
   /* mu = mu_c + mu_tilde*mu_tau; */
}

model{

  sigma ~ normal(0,10);

  /* mu_tilde ~ normal(0,1); */
  /* mu_tau ~ normal(0,1); */
  /* mu_c ~ normal(0,1); */
  alpha ~ normal(0,5);
  
  for(d in 1:D)
    {
      beta[d] ~ normal(0,5);
      //      alpha[d] ~normal(0,1);
    }

  #sigma ~normal(0,1);
  
  for(d in 1:D)
    {
      for(n in 2:N)
        X[n,d] ~ normal(X[n-1,d]*beta[d]+alpha,  sigma[d]);
    }

}

  

