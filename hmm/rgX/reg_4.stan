#include "../lib.stan"
#include "input.lstan"

transformed data{
  int K =3;
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
  real beta1;
  real beta2;
  real beta3;
  real<lower=0> kappa;
  simplex[K] rho;  
}

transformed parameters{
  vector[K] beta;

  beta[1]=beta1;
  beta[2]=beta2;
  beta[3]=beta3;
  
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

  beta1 ~ normal(1.5,3);
  beta3 ~ normal(0,2);
  beta2 ~ normal(-1.5,3);

  
  #include "reg_loglik.lstan"

}
#include "reg_gen.lstan"

