#include "../lib.stan"
#include "input.lstan"
parameters{
  
  //vector<lower=0> sigma;
  real<lower=0> sigma;
  simplex[K] tr[K];
  simplex[K] rho;
  vector[K] mu;
}

transformed parameters{

  matrix[K,K] Gamma= rep_matrix(0, K, K);
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }
}


model{
  sigma ~ normal(0,2);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho ~dirichlet(rep_vector(dir_rho/K,K));
  mu  ~normal(0,1);
  
  //  #include "normal_loglik.lstan"
}

//#include "normal_gen.lstan"

