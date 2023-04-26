#include "../lib.stan"
#include "input.lstan"

transformed data{
  int K =2;
}

parameters{
  //vector<lower=0>[K] sigma;

  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  //vector<lower=0>[K] sigma1;
  real<lower=0> sigma1;
  
  real<lower=0> mu1;
  real<upper=0> mu2;
  simplex[K] tr[K];
  simplex[K] rho;
  
}

transformed parameters{
  matrix[K,K] Gamma= rep_matrix(0, K, K);
  vector[K] mu;
  mu[1]=mu1;
  mu[2]=mu2;
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }
}


model{

  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));

  alpha0~normal(0,1);
  sigma1 ~normal(0,3);
  mu1 ~normal(0.5,1);
  mu2 ~normal(-0.5,1);

  #include "garch_loglik.lstan"
}

#include "garch_gen.lstan"
