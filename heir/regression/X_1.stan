data{

#include "../../input/X.stan"

}

transformed data{
  array[N] real x;
  for(n in 1:N)
    {
      x[n]=n*0.01;
    }
}

parameters{
  vector[D] beta;
  vector<lower=0>[D] sigma;
  //  vector[D] alpha;
  real alpha;
}

model{
  sigma ~ cauchy(0,5);
  alpha ~ normal(0,5);
  //  mu ~ normal(0,5);
  for(d in 1:D)
    {
      beta[d] ~ normal(0,5);
    }

  for(d in 1:D)
    {
      for(n in 1:N)
        {
          X[n,d]~ normal(x[n]*beta[d]+alpha,sigma[d]);
        }
    }
}

