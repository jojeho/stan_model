data{
  int T;
  int K;
  real y[T];
  real dir_tr;
  real dir_rho;
}

parameters{
  ordered [K] mu;
  vector[K] sigma;
  
  simplex[K] tr[K];
  simplex[K] rho;
}

transformed parameters{
  matrix[K,T] ob;
  matrix[K,K] Gamma= rep_matrix(0, K, K);

  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }

  
  for(k in 1:K)
    {
      for(t in 1:T)
        {
          ob[k,t]=normal_lpdf(y[t]|mu[k] ,sigma[k]);
        }
    }

}


model{
  
  sigma ~ normal(0,3);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  //rho  ~dirichlet(rep_vector(dir_rho/K,K));
  mu ~normal(0,2);

  target += hmm_marginal(ob,Gamma,rho);
}

generated quantities{
  matrix[K,T] prob = hmm_hidden_state_prob(ob,Gamma,rho);
  array[T] int label=hmm_latent_rng(ob,Gamma,rho);
  array[T] real y_hat;
  array[T] real log_lik;

  for( t in 1:T)
    {
      y_hat[t]=mu[label[t]];
      log_lik[t]=log_sum_exp(ob[,t]);
    }
}
