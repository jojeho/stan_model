
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
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  simplex[K] tr[K*N];
  simplex[K] rho;

}

transformed parameters{
  array[N] matrix[K,K] Gamma ;//= rep_matrix(0, K, K);
  for(n in 1:N)
    {
      for(k in 1:K)
        Gamma[n][k,]=tr[K*(n-1)+k]';
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
              ob[k,t]=normal_lpdf(yy[t]|mu[k] ,sigma);
            }
        }

      target += hmm_marginal(ob,Gamma[n],rho);
      pos +=T;
    }

}



