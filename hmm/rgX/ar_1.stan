#include "../lib.stan"
#include "input.lstan"

transformed data{
  int K =3;
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
  real alpha1;
  real alpha2;

}

transformed parameters{
  vector[K] beta;
  vector[K] alpha;
  beta[1]=beta1;
  beta[2]=0;
  beta[3]=beta2;
  alpha[1]=alpha1;
  alpha[2]=0;
  alpha[3]=alpha2;
  
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

  beta1 ~ von_mises(pi(),kappa); 
  beta2 ~ von_mises(-pi(),kappa);

  alpha1 ~ normal(0,1);
  alpha2 ~ normal(0,1);

  kappa~normal(5,2);

  #include "ar_loglik.lstan"

}
//#include "ar_gen.lstan"

