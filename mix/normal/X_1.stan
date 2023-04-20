data{
#include "../../input/X.stan"
  int K;
}


parameters{
  //ordered[K] alpha;
  ordered[K] mu;
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  simplex[K] theta;
}

transformed parameters{
  
}

model{

  sigma ~ cauchy(0,5);
  //alpha ~normal(0,1);
  mu ~normal(0,0.5);
  //  theta ~ dirichlet(rep_vector(0.3,K));

  
  for(d in 1:D)
    {
      vector[K] log_theta=log(theta);
      for( k in 1:K)
        {
          for(n in 1:N)
            {
              log_theta[k] +=normal_lpdf(X[n,d]|mu[k] ,sigma);
            }
        }
      target += log_sum_exp(log_theta);
    }
}


generated quantities{
  matrix[K,D] prob;
  vector[D] log_lik;
  
   for(d in 1:D)
    {
      vector[K] log_theta=theta;
      for( k in 1:K)
        {
          for(n in 1:N)
            {
              log_theta[k] *=exp(normal_lpdf(X[n,d]|mu[k] ,sigma));
            }
        }

      log_lik[d]=log_sum_exp(log_theta);
      for(k in 1:K)
        prob[k,d] =log_theta[k]/sum(log_theta);

    }
}


