#include "../lib.stan"
#include "input2.lstan"


/* transformed data{ */
/*   int K=6; */
/* } */

parameters{
  
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  simplex[K] tr[K];
  /* real<lower=0> mu1; */
  /* real<upper=0> mu2; */
  vector[K] mu;
  simplex[K] rho;
}

transformed parameters{
  matrix[K,K] Gamma= rep_matrix(0, K, K);
  
  
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }

}


model{
  sigma ~ normal(0,1);
  mu ~ normal(0,0.3);
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
    
#include "normal_loglik.lstan"

}

#include "normal_gen.lstan"


