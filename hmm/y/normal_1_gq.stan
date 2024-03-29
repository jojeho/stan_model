data{
  int T;
  int K;
  real y[T];
}

parameters{
  ordered [K] mu;
  vector<lower=0>[K] sigma;
  
  simplex[K] tr[K];
  simplex[K] rho;
}

/* model{ */
/*   sigma ~ normal(0,0.1); */
/*   tr  ~ dirichlet(rep_vector(dir_tr/(K),K)); */
/*   rho  ~dirichlet(rep_vector(dir_rho/K,K)); */
/*   mu ~normal(0,0.01); */

/* } */

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


generated quantities{
  matrix[K,T] prob = hmm_hidden_state_prob(ob,Gamma,rho);

  vector[K] n_prob=rep_vector(0,K);

  for(k in 1:K)
    {
      for(j in 1:K)
        {
          n_prob[k] +=log(prob[j,T])+log(Gamma[j,k]);
        }
    }
  n_prob=softmax(n_prob);
  
}
