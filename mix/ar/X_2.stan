data{
#include "../../input/X.stan"
  vector[N] x;
  int K;
}


parameters{
  real alpha;
  ordered[K] beta;
  vector<lower=0>[K] sigma;
  //real<lower=0> sigma;
  simplex[K] theta;
}

transformed parameters{
  
}

model{

  sigma ~ cauchy(0,5);
  alpha ~normal(0,1);
  beta ~normal(0,5);
  //theta ~ dirichlet(rep_vector(10/K,K));

  for(d in 1:D)
    {
      vector[K] log_theta=log(theta);
      for( k in 1:K)
        {
          for(n in 1:N)
            {
              log_theta[k] +=normal_lpdf(X[n,d]|exp(beta[k]*x[n]+alpha) ,sigma[k]);
            }
        }
      target += log_sum_exp(log_theta);
    }
}





/* generated quantities{ */
/*   matrix[K,D] prob; */
    
/*    for(d in 1:D) */
/*     { */
/*       vector[K] log_theta=theta; */
/*       for( k in 1:K) */
/*         { */
/*           for(n in 1:N) */
/*             { */
/*               log_theta[k] *=exp(normal_lpdf(X[n,d]|exp(beta[k]*x[n]) ,sigma)); */
/*             } */
/*         } */
/*       for(k in 1:K) */
/*         prob[k,d] =log_theta[k]/sum(log_theta); */

/*     } */
/* } */


