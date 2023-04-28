#include "../lib.stan"
#include "input.lstan"


transformed data{
  int K=3;
}

parameters{
  
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  simplex[K] tr[K];
  /* real<lower=0> mu1; */
  /* real<upper=0> mu2; */
  real mu1;
  real mu2;
  real mu3;
  
  simplex[K] rho;
}

transformed parameters{
  matrix[K,K] Gamma= rep_matrix(0, K, K);
  vector[K] mu;
  
  mu[1]=mu1;
  mu[2]=mu2;
  mu[3]=mu3;
  
  
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }

}


model{
  sigma ~ normal(0,0.2);
  mu1 ~normal(0.2,0.1);
  mu2 ~normal(0,0.1);
  mu3 ~normal(-0.2,0.1);
  rho  ~dirichlet(rep_vector(dir_rho/K,K));

#include "normal_loglik.lstan"

}

#include "normal_gen.lstan"


