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

  vector[D] beta;
  vector[D] beta1;
}

transformed parameters{

}

model{


  alpha ~ normal(0,1);
  //  tau ~ cauchy(0,5);
  beta~normal(0,2);
  beta1~normal(0,2);

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
