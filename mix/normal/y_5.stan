data{
#include "../../input/y.stan"
  int K;
}


parameters{
  //  ordered[K] alpha;
  ordered[K] mu;
  real<lower=0> sigma;
  simplex[K] theta;
}


model{

  sigma ~ cauchy(0,5);
  //  alpha ~normal(0,0.1);
  mu ~normal(0,5);
  // theta ~ dirichlet(rep_vector(2,K));

  matrix[K,T] ob;
  for(t in 1:T)
    {
      vector[K] log_theta=log(theta);
      for( k in 1:K)
        log_theta[k] +=normal_lpdf(y[t]|mu[k] ,sigma);
      target += log_sum_exp(log_theta);
    }
}


generated quantities{
  array[T] vector[K] prob;
  array[T] real y_hat;
  array[T] real log_lik;
  for(t in 1:T)
    {
      vector[K] log_theta=log(theta);
      for( k in 1:K)
        log_theta[k] +=normal_lpdf(y[t]|mu[k] ,sigma[k]);
      prob[t]=softmax(log_theta);
      y_hat[t]=normal_rng(dot_product(prob[t],to_vector(mu)),dot_product(prob[t],sigma));
      log_lik[t]=normal_lpdf(y[t]|dot_product(prob[t],to_vector(mu)),dot_product(prob[t],sigma)));
    }
}

