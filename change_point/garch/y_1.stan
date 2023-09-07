data{
  #include "../../input/y.stan"
}


transformed data {
  real log_unif;
  real alpha_prior;
  real sigma_prior;
  real mu_prior;
  log_unif = -log(T);
  mu_prior=0.1;
  alpha_prior=0.1;
  sigma_prior=1;
    
}
parameters {
  real mu1;
  real mu2;

  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;

  real<lower=0> sigma1_1;
  real<lower=0> sigma1_2;
}

transformed parameters {

  vector[T] lp;
  vector<lower=0>[T] sigma_1;
  vector<lower=0>[T] sigma_2;
  
  {
    sigma_1[1] = sigma1_1;
    sigma_2[1] = sigma1_2;
    
    for (t in 2:T)
      {
        sigma_1[t] = sqrt(alpha0
                          + alpha1 * pow(y[t-1] - mu1, 2)
                          + beta1 * pow(sigma_1[t-1], 2));

        sigma_2[t] = sqrt(alpha0
                          + alpha1 * pow(y[t-1] - mu2, 2)
                          + beta1 * pow(sigma_2[t-1], 2));        
      }

    
    vector[T + 1] lp_e;
    vector[T + 1] lp_l;
    lp_e[1] = 0;
    lp_l[1] = 0;

    for (t in 1:T) {
      lp_e[t + 1] = lp_e[t] + normal_lpdf(y[t] | mu1,sigma_1[t]);
      lp_l[t + 1] = lp_l[t] + normal_lpdf(y[t] | mu2,sigma_2[t]);
    }
    lp = rep_vector(log_unif + lp_l[T + 1], T)
      + head(lp_e, T) - head(lp_l, T);
  }   
}
model {
  
  alpha1 ~ normal(0,1);

  alpha0 ~ normal(0,alpha_prior);
  
  beta1 ~ normal(0,1);
  
  sigma1_1 ~ normal(0,sigma_prior);
  sigma1_2 ~ normal(0,sigma_prior);
  mu1 ~ normal(0,mu_prior);
  mu2 ~ normal(0,mu_prior);

  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  array[T] real log_lik;
  array[T] real y_hat;
  s = categorical_logit_rng(lp);
  for(t in 1:T)
    {
      if( t < s)
        {
          log_lik[t]=normal_lpdf(y[t]|mu1,sigma_1[t]);
          y_hat[t]=normal_rng(mu1,sigma_1[t]);
        }
      else
        {
          log_lik[t]=normal_lpdf(y[t]|mu2,sigma_2[t]);
          y_hat[t]=normal_rng(mu2,sigma_2[t]);
        }
    }
  
}
