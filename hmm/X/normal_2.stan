#include "../lib.stan"

data{
  int N;
  int D;
  array[N] vector[D] X;
  real dir_tr;
  real dir_rho;
}

transformed data{
  int K=5;
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
  for(k in 1:K)
    {
      rho[k]=0.2;
    }
  
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }

  matrix[K,D] ob=rep_matrix(0,K,D);
  for(d in 1:D)
    {
      for(k in 1:K)
        {
          for(n in 1:N)
            {
              ob[k,d] +=normal_lpdf(X[n,d]|mu[k] ,sigma);
            }
        }
    }
}


model{
  sigma ~ normal(0,1);
  mu1 ~student_t(5,-4,2);
  mu2 ~normal(-2,0.5);
  mu3 ~normal(0,0.5);
  mu4 ~normal(2,0.5);
  mu5 ~student_t(5,4,2);
  target += hmm_marginal(ob,Gamma,rho);
}

generated quantities{
  //matrix[K,D] prob=hmm_hidden_state_prob(ob,Gamma,rho);
  matrix[K,D] prob=hmm_forward_prob(ob,Gamma,rho);
}
