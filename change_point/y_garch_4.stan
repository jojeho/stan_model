data{
  #include "../input/y.stan"
}


transformed data {
  real log_unif;
  real alpha_prior;
  real sigma_prior;
  real mu_prior;
  log_unif = -log(T);
  mu_prior=0.3;
  alpha_prior=0.3;;
  sigma_prior=1;
}

parameters {
  real mu;

  real<lower=0> alpha0_1;
  real<lower=0,upper=1> alpha1_1;
  real<lower=0,upper=(1-alpha1_1)> beta1_1;
  
  real<lower=0> alpha0_2;
  real<lower=0,upper=1> alpha1_2;
  real<lower=0,upper=(1-alpha1_2)> beta1_2;
  
  real<lower=0> sigma1_1;
  real<lower=0> sigma1_2;

  real mu_beta;
}

transformed parameters {
  real mu1;
  real mu2;
  vector[T] lp;
  vector<lower=0>[T] sigma_1;
  vector<lower=0>[T] sigma_2;

  mu1=mu;
  mu2=mu1*mu_beta;
  
  {
    sigma_1[1] = sigma1_1;
    sigma_2[1] = sigma1_2;
    
    for (t in 2:T)
      {
        sigma_1[t] = sqrt(alpha0_1
                          + alpha1_1 * pow(y[t-1] - mu1, 2)
                          + beta1_1 * pow(sigma_1[t-1], 2));

        sigma_2[t] = sqrt(alpha0_2
                          + alpha1_2 * pow(y[t-1] - mu2, 2)
                          + beta1_2 * pow(sigma_2[t-1], 2));        
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
  
  alpha1_1 ~ normal(0,1);
  alpha1_2 ~ normal(0,1);

  alpha0_1 ~ normal(0,alpha_prior);
  alpha0_2 ~ normal(0,alpha_prior);
  
  beta1_1 ~ normal(0,1);
  beta1_2 ~ normal(0,1);
  
  sigma1_1 ~ normal(0,sigma_prior);
  sigma1_2 ~ normal(0,sigma_prior);
  mu ~ normal(0,mu_prior);
  mu_beta ~ normal(-1,0.5);

  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  array[T] real log_lik;  
  s = categorical_logit_rng(lp);
}
