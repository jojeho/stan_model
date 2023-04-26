#include "../lib.stan"
#include "input.lstan"

transformed data{
  int K =2;
  vector[N] x;
  for(t in 1:N)
    {
      x[t]=t*0.1;
    }  
}

parameters{
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  simplex[K] tr[K];
  //  simplex[K] rho;
  real<lower=0> beta1;
  real<upper=0> beta2;
  real<lower=0> kappa;
  simplex[K] rho;  
}

transformed parameters{
  vector[K] beta;

  beta[1]=beta1;
  beta[2]=beta2;
  
  matrix[K,K] Gamma= rep_matrix(0, K, K);
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }
}


model{
  sigma ~ normal(0,1);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  beta1 ~ von_mises(pi()*2,kappa); 
  beta2 ~ von_mises(-pi()*2,kappa);
  kappa~normal(5,2);

  #include "reg_loglik.lstan"
}
#include "reg_gen.lstan"
