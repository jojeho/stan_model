data{
#include  "input.lstan"
  int K;
}


parameters{
  real alpha;
  real beta;
  //vector<lower=0>[K] sigma;
  positive_ordered[K] sigma;
  //real<lower=0> sigma1;
  ///  real<lower=0> sigma2;
  simplex[K] theta;
}

transformed parameters{
  /* vector[K] sigma; */
  /* sigma[1]=sigma1; */
  /* sigma[2]=sigma2; */
}

model{

  sigma ~ cauchy(0.1,2);
  //  sigma2 ~ cauchy(0.3,1);
  //alpha ~normal(0,1);
  beta ~normal(0,2);
  alpha ~normal(0,1);
  //theta ~ dirichlet(rep_vector(0.3,K));

  
  for(d in 1:D)
    {
      vector[K] log_theta=log(theta);
      for( k in 1:K)
        {
          for(n in 2:N)
            {
              log_theta[k] +=normal_lpdf(X[n,d]|X[n-1,d]*beta+alpha,sigma[k]);
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
          for(n in 2:N)
            {
              log_theta[k] *=exp(normal_lpdf(X[n,d]|X[n-1,d]*beta+alpha ,sigma[k]));
            }
        }

      log_lik[d]=log_sum_exp(log_theta);
      for(k in 1:K)
        prob[k,d] =log_theta[k]/sum(log_theta);

    }
}


