
data{
#include "../../input/X.stan"
  vector[N] x;
  int K;
}


parameters{
  vector[K] alpha;
  ordered[K] beta;
  vector<lower=0>[K] sigma;
  simplex[K] theta;
}


model{

  sigma ~ cauchy(0,5);
  alpha ~normal(0,1);
  beta ~normal(0,1);
  //theta ~ dirichlet(rep_vector(3/K,K));

  for(d in 1:D)
    {
      vector[K] log_theta=log(theta);
      for( k in 1:K)
        {
          for(n in 1:N)
            {
              log_theta[k] +=normal_lpdf(X[n,d]|alpha[k]+beta[k]*x[n] ,sigma[k]);
            }
        }
      target += log_sum_exp(log_theta);
    }
}


/* generated quantities{ */
/*   matrix[K,D] prob; */
    
/*    for(d in 1:D) */
/*     { */
/*       vector[K] log_theta=log(theta); */
/*       for( k in 1:K) */
/*         { */
/*           for(n in 1:N) */
/*             { */
/*               log_theta[k] +=normal_lpdf(X[n,d]|alpha[k]+beta[k]*x[n] ,sigma[k]); */
/*             } */
/*         } */
/*       prob[,d] =softmax(log_theta); */
/*     } */
/* }  */

