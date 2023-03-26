data{

#include "../../input/X.stan"
  array[N] real x;
}

parameters{
  vector[D] beta;
  vector<lower=0>[D] sigma;
  vector[D] alpha;
}

model{

  sigma ~ cauchy(0,5);
  alpha ~ normal(0,5);
  beta  ~ normal(0,5);

  for(d in 1:D)
    {
      for(n in 1:N)
        {
          X[n,d]~ normal(x[n]*beta[d]+alpha[d],sigma[d]);
        }
    }
}

