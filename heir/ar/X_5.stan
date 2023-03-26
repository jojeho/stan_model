data{

#include "../../input/X.stan"

  array[N] real x;
}


parameters{
  
  vector<lower=0>[D] sigma;
  //  real alpha;

  //vector[D] beta1;
  //real alpha_tilde;
  //real tau;
  real alpha;

  real beta_tau;
  vector[D] beta_tilde;
  real beta_c;

  real beta1_tau;
  vector[D] beta1_tilde;
  real beta1_c;
}

transformed parameters{
  vector[D] beta;
  vector[D] beta1;
  for(d in 1:D)
    {
      beta[d]=beta_c + beta_tilde[d]*beta_tau;
      beta1[d]=beta1_c + beta1_tilde[d]*beta1_tau;
    }
}

model{

  beta_tilde ~ normal(0,1);
  beta_tau~cauchy(0,2);
  beta_c ~ normal(0,2);

  beta1_tilde ~ normal(0,1);
  beta1_tau~cauchy(0,2);
  beta1_c ~ normal(0,2);


  alpha ~ normal(0,1);
  //  tau ~ cauchy(0,5);
  //  beta~normal(0,2);
  //beta1~normal(0,2);

  for(d in 1:D)
    {
      array[N] real y_hat;
      for(n in 1:N)
        {
          //y_hat[t] = beta1*exp(-x[t]*beta+gamma)+alpha;
          y_hat[n] = beta1[d]*exp(-x[n]*beta[d])+alpha;
        }
      for(n in 1:N)
        {
          X[n,d] ~ normal(y_hat[n],sigma[d]);
        }
    }
  
}
