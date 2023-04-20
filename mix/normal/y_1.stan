data{
#include "../../input/y.stan"
  int K;
}


parameters{
  //  ordered[K] alpha;
  ordered[K] mu;
  vector<lower=0>[K] sigma;
  simplex[K] theta;
}


model{

  sigma ~ cauchy(0,0.5);
  //  alpha ~normal(0,0.1);
  mu ~normal(0,0.1);
  //theta ~ dirichlet(rep_vector(2,K));

    int pos=1;
    matrix[K,T] ob;
    for(t in 1:T)
      {
          vector[K] log_theta=log(theta);
          for( k in 1:K)
            log_theta[k] +=normal_lpdf(y[t]|mu[k] ,sigma[k]);
          target += log_sum_exp(log_theta);
      }
}


/* generated quantities{ */
/*   array[T] vector[K] prob; */
/*   for(t in 1:T) */
/*     { */
/*       vector[K] log_theta=log(theta); */
/*       for( k in 1:K) */
/*         log_theta[k] +=normal_lpdf(y[t]|mu[k] ,sigma[k]); */
/*       prob[t]=softmax(log_theta); */
/*     } */
/* } */

