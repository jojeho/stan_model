#include "../lib.stan"

data{
  int N;
  int D;
  array[N] vector[D] X;
  real dir_tr;
  real dir_rho;
  //  real sigma_prior;    
  int K;
}


parameters{
  ordered[K] mu;
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
  
  matrix[K,K] Gamma= rep_matrix(0, K, K);
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }
  matrix[K,D] ob=rep_matrix(0,K,D);
  for(d in 1:D)
    {
      for(n in  1:N)
        {
          for(k in 1:K)
            ob[k,d] =normal_lpdf(X[n,d]|mu[k] ,sigma);
        }
    }
  
}


model{
  sigma ~ normal(0,1);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  beta1 ~ von_mises(pi()*2,kappa); 
  beta2 ~ von_mises(-pi()*2,kappa);
  kappa~normal(5,2);
  target += hmm_marginal(ob,Gamma,rho);
}

generated quantities{
  matrix[K,D] prob=hmm_forward_prob(ob,Gamma,rho);
}
