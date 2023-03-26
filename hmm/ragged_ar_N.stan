
data{
  int N;
  int AT;
  int K;
  array[AT] real  y;
  array[N] int GT;
  array[N] real py;
  
  real dir_tr;
  real dir_rho;

}

parameters{
  ordered[K] beta;
  vector[K] alpha;
  vector<lower=0>[K] sigma;
  simplex[K] tr[K];
  simplex[K] rho;

}

transformed parameters{
  matrix[K,K] Gamma= rep_matrix(0, K, K);
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }
}


model{

  sigma ~ normal(0,0.1);
  
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
#  for( n in 1:N)
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  
  alpha ~normal(0,0.1);
  beta ~normal(1,0.1);

  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];

      matrix[K,T] ob;

      array[T] real yy=segment(y,pos,T);

      for(k in 1:K)
        {
          ob[k,1]=normal_lpdf(yy[1]|alpha[k]+beta[k]*py[n] ,sigma[k]);
          for(t in 2:T)
            {
              ob[k,t]=normal_lpdf(yy[t]|alpha[k]+beta[k]*yy[t-1] ,sigma[k]);
            }
        }

      target += hmm_marginal(ob,Gamma,rho);
      pos +=T;
    }
  


}

/* generated quantities{ */
/*   vector[AT] oblik; */
/*   matrix[K,AT] prob; */

/*   int pos=1; */
/*   for(n in 1:N) */
/*     { */
/*       int T=GT[n]; */
/*       matrix[K,T] ob; */
/*       matrix[K,T] local_prob; */
/*       array[T] real yy=segment(y,pos,T); */
  
/*       for(k in 1:K) */
/*         { */
/*           for(t in 1:T) */
/*             { */
/*               ob[k,t]=normal_lpdf(yy[t]|mu[k] ,sigma[k]); */
/*             } */
/*         } */

/*       local_prob=hmm_hidden_state_prob(ob,Gamma,rho[n]); */
      
/*       for( t in 1:T) */
/*         { */
/*           oblik[pos+t-1]=log_sum_exp(ob[,t]); */
/*           prob[:,pos+t-1]=local_prob[:,t]; */
/*         } */
      
/*       pos +=T; */
/*     } */
  
/* } */
