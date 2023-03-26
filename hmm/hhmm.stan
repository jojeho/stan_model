data{
  int T;
  array[T] int hy;
  array[T] real x;
  int hK ;
  int K;
  real dir_tr;
  real dir_rho;
  int i_hy;
}


parameters{
  ordered [K] mu;
  //vector<lower=0>[hK] sigma;
  real<lower=0> sigma;
  simplex[K] tr[K];
  simplex[K] rho;
}


transformed parameters{
  matrix[K,K] Gamma= rep_matrix(0, K, K);
  array[T] real ta;
  matrix[K,T] c_oblik;
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }
  ta[1] =log(Gamma[i_hy,hy[1]]);
  for(t in 2:T)
    {
      ta[t]=log(Gamma[hy[t-1],hy[t]]);
    }
  for(t in 1:T)
    {
      for(k in 1:K)
        {
          c_oblik[k,t]=normal_lpdf(x[t]|c_mu[k],c_sigma);
        }
    }

  for (t in 2:T) {
    for (j in 1:K) { // j = current (t)
      for (i in 1:K) { // i = previous (t-1)
            // Murphy (2012) Eq. 17.48
            // belief state + transition prob + local evidence at t
            accumulator[i] = unalpha_tk[t-1, i] + log(A[j,i]) + c_oblik[j,t]
              }

      unalpha_tk[t, j] = log_sum_exp(accumulator);
    }
  }

  for(t in 2:T)
    {
      for(j in 1:K)
        {
          for(i in 1:K)
            {
              accumulator[i] = hunalpha_tk[t, i-1] + log(MA[j,i]) + H_oblik[j,t] 
            }

          unalpha_tk[t, j] = log_sum_exp(accumulator);
        }

      for(j in 1:K)
        {
          for(i in 1:K)
            {
              accumulator[i] = log(HMA[j,i])+unalpha_tk[t,i];
            }
          hunalpha_tk[t,j]=log_sum_exp(accumulator);
        }
    }

}

model{
  sigma ~ normal(0,3);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  mu ~normal(0,10);
  for(t in 1:T)
    {
      target+= ta[t];
      target += log_sum_exp(hunalpha_tk[T]); // Note: update based only on last unalpha_tk
    }
  
}

generated quantities{
  matrix[K,T] oblik;
  for(t in 1:T)
    {
      for(k in 1:K)
        {
          oblik[k,t]=normal_lpdf(x[t]|mu[k],sigma);
        }
    }
  
  matrix[K,T] prob=hmm_hidden_state_prob(oblik,Gamma,rho); 
 } 
