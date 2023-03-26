
data{
  int N;
  int AT;
  array[AT] real  y;
  array[N] int GT;
  int K;
  real dir_tr;
  real dir_rho;
}

parameters{
  ordered[K] mu;
  vector<lower=0>[K] sigma;
  //real<lower=0> sigma;
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
  sigma ~ cauchy(0,3);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  mu ~normal(0,1);

  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];

      matrix[K,T] ob;

      array[T] real yy=segment(y,pos,T);
      
  
      for(k in 1:K)
        {
          for(t in 1:T)
            {
              ob[k,t]=normal_lpdf(yy[t]|mu[k] ,sigma[k]);
            }
        }

      target += hmm_marginal(ob,Gamma,rho);
      pos +=T;
    }

}

/* generated quantities{ */
/*   matrix[K,NT] prob; */
/*   //array[N] matrix[K,77] prob; */
  
/*   int pos=1; */
/*   for(n in 1:NN) */
/*     { */
/*       int T=NGT[n]; */

/*       matrix[K,T] n_ob; */
/*       matrix[K,T] local_prob; */
/*       array[T] real ny=segment(n_y , pos ,T); */

/*       for(t in 1:T) */
/*         { */
/*           for(k in 1:K) */
/*             { */
/*               n_ob[k,t]=normal_lpdf(ny[t]|mu[k],sigma[k]); */
/*             } */
/*         } */
/*       local_prob = hmm_hidden_state_prob(n_ob,Gamma,rho); */
/*       for( t in 1:T) */
/*         { */
/*           prob[:,pos+t-1]=local_prob[:,t]; */
/*           // prob[n][:,t]=local_prob[:,t]; */
/*         } */
/*       pos +=T; */
/*     } */

/* } */
