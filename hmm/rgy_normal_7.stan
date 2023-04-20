
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
  //  real<lower=0> sigma;
  simplex[K] tr[K*N];
  simplex[K] rho;
  real<lower=0> sigma_beta;
      
}

transformed parameters{
  vector<lower=0>[K] sigma;
  array[N] matrix[K,K] Gamma ;//= rep_matrix(0, K, K);
  for(n in 1:N)
    {
      for(k in 1:K)
        Gamma[n][k,]=tr[K*(n-1)+k]';
    }
  for(k in 1:K)
    {
      sigma[k]=mu[k]^2*sigma_beta;
    }

}

model{
  //sigma ~ cauchy(0,3);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  mu ~normal(0,3);
  sigma_beta ~ normal(0,2);

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

      target += hmm_marginal(ob,Gamma[n],rho);
      pos +=T;
    }

}



