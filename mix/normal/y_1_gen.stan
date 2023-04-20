data{
#include "../../input/y.stan"
  int K;
  vector[K] mu_s;
  vector[K] sigma;
  vector[K] theta_s;
}

transformed data{
  vector[K] mu;
  vector[K] theta;
  mu=mu_s/sum(mu_s);
  theta=theta_s/sum(theta_s);
    
}

generated quantities{
  array[T] vector[K] prob;
  for(t in 1:T)
    {
      vector[K] log_theta=log(theta);
      for( k in 1:K)
        log_theta[k] +=normal_lpdf(y[t]|mu[k] ,sigma);
      prob[t]=softmax(log_theta);
    }
}

