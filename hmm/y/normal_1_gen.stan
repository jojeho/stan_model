data{
  int T;
  int K;
  real y[T];
  real dir_tr;
  real dir_rho;
  real mu_scale_prior;
}

parameters{
  ordered [K] mu;
  vector<lower=0>[K] sigma;
  
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

  //  sigma ~ gamma(1,0.1);
  sigma ~ normal(0,10);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  mu ~normal(0,mu_scale_prior);
  

  target += hmm_marginal(ob,Gamma,rho);
}


generated quantities{
  matrix[K,T] prob = hmm_hidden_state_prob(ob,Gamma,rho);
  vector[K] n_prob=rep_vector(0,K);
  vector[T] oblik;

  for(k in 1:K)
    {
      for(j in 1:K)
        {
          n_prob[k] +=log(prob[j,T])+log(Gamma[j,k]);
        }
    }
  n_prob=softmax(n_prob);
  
  for(t in 1:T)
    {
      oblik[t]=log_sum_exp(log(prob[,t])+ob[,t]);
    }
}
