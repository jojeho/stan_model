
data{
  int N;
  int AT;
  int K;
  array[AT] real  y;
  array[N] int GT;

}


parameters{
  ordered[K] alpha;
  real beta;
  vector<lower=0>[K] sigma;
  //real<lower=0> sigma;
  simplex[K] theta;
}

transformed parameters{
  
}

model{

  sigma ~ normal(0,3);
  alpha ~normal(0,2);
  beta ~normal(0,5);
  #theta ~ dirichlet(rep_vector(2/K,K));

  
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];

      matrix[K,T] ob;

      array[T] real yy=segment(y,pos,T);

      for(t in 2:T)
        {
          vector[K] log_theta=log(theta);
          for( k in 1:K)
            log_theta[k] +=normal_lpdf(yy[t]|alpha[k]+beta*yy[t-1] ,sigma[k]);
          target += log_sum_exp(log_theta);
        }
      
      pos +=T;
    }
  
}

