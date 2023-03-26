data{
  int T;
  int K;
  real y[T];
  real dir_tr;
  real dir_rho;
}

parameters{
  ordered [K] mu;
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  
  simplex[K] tr[K];
  simplex[K] rho;
}

transformed parameters{
  matrix[K,T] ob;
    
  for(k in 1:K)
    {
      for(t in 1:T)
        {
          ob[k,t]=normal_lpdf(y[t]|mu[k] ,sigma);
        }
    }

  

  matrix[K,K] Gamma= rep_matrix(0, K, K);

  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }

  }


model{
  
  sigma ~ cauchy(0,3);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  mu ~normal(0,2);

  target += hmm_marginal(ob,Gamma,rho);
}

/* generated quantities{ */
/*   //matrix[K,T] prob=hmm_hidden_state_prob(ob,Gamma,rho); */
/*   array[T] int label=hmm_latent_rng(ob,Gamma,rho); */
/*   vector[T] y_hat; */
/*   vector[T] log_lik; */
/*   for(t in 1:T) */
/*     { */
/*       y_hat[t]=normal_rng(mu[label[t]],sigma); */
/*       log_lik[t]=normal_lpdf(y[t]|mu[label[t]],sigma); */
/*     } */
/* } */



