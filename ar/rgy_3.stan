data{

#include "../input/rgy.stan"
  real PM_alpha;
  real PM_beta;
  real PM_sigma;

  real PS_alpha;
  real PS_beta;
  real PS_sigma;


}

parameters{
  array[N] real alpha;
  array[N] real beta;
  real<lower=0> sigma;
}

model{
  for(n in 1:N)
    {
      beta[n] ~ normal(PM_beta,PS_beta);
      alpha[n] ~normal(PM_alpha,PS_alpha);
    
    }

  sigma ~normal(PM_sigma,PS_sigma);
  
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      array[T] real yy=segment(y,pos,T);
      for(t in 2:T)
        {
          yy[t] ~ normal(yy[t-1]*beta[n]+alpha[n] , sigma);
        }
      pos +=T;
    }
}

generated quantities{
  array[AT] real y_hat;//=rep_array(0,AT);
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      array[T] real yy=segment(y,pos,T);
      for(t in 2:T)
        y_hat[pos+t-1] = normal_rng(yy[t-1]*beta[n]+alpha[n] , sigma);
      pos +=T;
    }
}
  

