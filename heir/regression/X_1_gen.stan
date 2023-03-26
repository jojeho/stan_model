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

generated quantities{
  array[N] vector[D] X_hat;
  
  for(d in 1:D)
    {
      for(n in 1:N)
        {
          X_hat[n,d]= normal_rng(x[n]*beta[d]+alpha,sigma[d]);
        }
    }
}
