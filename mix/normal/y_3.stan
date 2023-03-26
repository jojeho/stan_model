data{
#include "../../input/y.stan"
  int K;
}


parameters{
  //  ordered[K] alpha;
  real mu;
  positive_ordered[K] sigma;
  simplex[K] theta;
}


model{

  sigma ~ normal(0,0.2);
  //  alpha ~normal(0,0.1);
  mu ~normal(1,0.1);
  //  theta ~ dirichlet(rep_vector(2/K,K));

    int pos=1;
    matrix[K,T] ob;
    for(t in 1:T)
      {
          vector[K] log_theta=log(theta);
          for( k in 1:K)
            log_theta[k] +=normal_lpdf(y[t]|mu ,sigma[K]);
          target += log_sum_exp(log_theta);
      }
}

generated quantities{
  array[T] vector[K] prob;
  for(t in 1:T)
    {
      vector[K] log_theta=log(theta);
      for( k in 1:K)
        log_theta[k] +=normal_lpdf(y[t]|mu ,sigma[k]);
      prob[t]=softmax(log_theta);
    }
}
