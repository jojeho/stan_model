data{
  #include "../input/rgy.stan"
}

transformed data {
  real alpha_prior;
  real sigma_prior;
  real mu_prior;
  mu_prior=0.1;
  alpha_prior=0.1;
  sigma_prior=1;
    
}


parameters {
  array[N] real  mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  //array[N] real<lower=0> sigma1;
  vector<lower=0>[N] sigma1;
  
}

model {

  mu ~ normal(0,mu_prior);
  alpha1 ~ normal(0,1);
  alpha0 ~ normal(0,alpha_prior);
  beta1 ~ normal(0,1);
  sigma1 ~ normal(0,sigma_prior);

  
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      vector[T] yy=segment(y,pos,T);
      array[T] real sigma;
      sigma[1]=sigma1[n];
      for (t in 2:T)
        {
          sigma[t] = sqrt(alpha0
                          + alpha1 * pow(yy[t-1] - mu[n], 2)
                          + beta1 * pow(sigma[t-1], 2));
        }
      
      yy~ normal(mu[n], sigma);      
      
      pos +=T;
    }
  
}

