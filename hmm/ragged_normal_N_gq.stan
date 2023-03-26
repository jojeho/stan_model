
data{
  int N;
  int AT;
  int K;
  array[AT] real  y;
  array[N] int GT;
  
  real dir_tr;
  real dir_rho;

}

parameters{
  ordered[K] mu;
  vector<lower=0>[K] sigma;
  simplex[K] tr[K];
  array[N] simplex[K] rho;

}

transformed parameters{
  matrix[K,K] Gamma= rep_matrix(0, K, K);
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }
}



generated quantities{
  vector[AT] log_lik;
  matrix[K,AT] prob;

  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      matrix[K,T] ob;
      matrix[K,T] local_prob;
      array[T] real yy=segment(y,pos,T);
  
      for(k in 1:K)
        {
          for(t in 1:T)
            {
              ob[k,t]=normal_lpdf(yy[t]|mu[k] ,sigma[k]);
            }
        }

      local_prob=hmm_hidden_state_prob(ob,Gamma,rho[n]);
      
      for( t in 1:T)
        {
          log_lik[pos+t-1]=log_sum_exp(ob[,t]);
          prob[:,pos+t-1]=local_prob[:,t];
        }
      
      pos +=T;
    }
  
}
