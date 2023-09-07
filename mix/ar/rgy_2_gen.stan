
data{
  int N;
  int AT;
  int K;
  array[AT] real  y;
  array[N] int GT;

}


parameters{
  //  ordered[K] alpha;
  ordered[K] beta;
  vector<lower=0>[K] sigma;
  simplex[K] theta;
}


model{

  sigma ~ normal(0,0.1);
  //  alpha ~normal(0,0.1);
  beta ~normal(1,0.1);
  #theta ~ dirichlet(rep_vector(2/K,K));

  
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];

      matrix[K,T] ob;

      array[T] real yy=segment(y,pos,T);
      vector[K] log_theta=log(theta);

      for(t in 2:T)
        {
          for( k in 1:K)
            log_theta[k] +=normal_lpdf(yy[t]|beta[k]*yy[t-1] ,sigma[k]);
        }
      target += log_sum_exp(log_theta);      
      pos +=T;
    }
  
}

generated quantities{

  array[N] vector[K] prob;
  {
    int pos=1;
    for(n in 1:N)
      {
        int T=GT[n];

        matrix[K,T] ob;

        array[T] real yy=segment(y,pos,T);
        vector[K] log_theta=log(theta);
        for(t in 2:T)
          {
            for( k in 1:K)
              {
                log_theta[k] +=normal_lpdf(yy[t]|beta[k]*yy[t-1] ,sigma[k]);
              }


        prob[n]=softmax(log_theta);
        pos +=T;
      }
  }
}
