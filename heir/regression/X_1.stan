data{

#include "../../input/X.stan"

}

transformed data{
  array[N] real x;
  for(n in 1:N)
    {
      x[n]=n*0.1;
    }
}

parameters{
  vector[D] beta;
  vector<lower=0>[D] sigma;
}

model{
  sigma ~ cauchy(0,5);
  //  mu ~ normal(0,5);
  for(d in 1:D)
    {
      beta[d] ~ normal(0,5);
    }

  for(d in 1:D)
    {
      for(n in 1:N)
        {
          X[n,d]~ normal(x[n]*beta[d],sigma[d]);
        }
    }
}

