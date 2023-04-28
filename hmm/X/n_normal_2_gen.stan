#include "../lib.stan"
#include "n_normal_gen_input.stan"


transformed data{
  int K=7;
}

parameters{
  
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  simplex[K] tr[K];
  real mu1;
  real mu2;
  real mu3;
  real mu4;
  real mu5;
  real mu6;
  real mu7;
}

transformed parameters{
  matrix[K,K] Gamma= rep_matrix(0, K, K);
  vector[K] mu;
  simplex[K] rho;
  mu[1]=mu1;
  mu[2]=mu2;
  mu[3]=mu3;
  mu[4]=mu4;
  mu[5]=mu5;
  mu[6]=mu6;
  mu[7]=mu7;
  
  rho[1]=0.05;
  rho[2]=0.1;
  rho[3]=0.2;
  rho[4]=0.3;
  rho[5]=0.2;
  rho[6]=0.1;
  rho[7]=0.05;
  
  
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }

}


#include "n_normal_gen_lib.stan"

/* generated quantities{ */
/*   //matrix[K,D] prob=hmm_hidden_state_prob(ob,Gamma,rho); */
/*   matrix[K,D] prob=hmm_forward_prob(ob,Gamma,rho); */
/* } */

