data{
 #include "../../input/X.stan"
}

parameters{
  vector[D] mu;
  real mu_alpha;
  vector<lower=0>[D] sigma;
  real<lower=0> ps_mu;
  real<lower=0> ps_sigma;
  //real<lower=0> sigma_alpha;
}

model{

  #ps_mu~normal(0,1);
  mu ~cauchy(0,10);
  ps_sigma ~ inv_gamma(0.3,0.3);
  ps_mu ~ normal(0,10);
  mu_alpha ~normal(0,10);
  //  sigma_alpha ~cauchy(0,10);
  for(d in 1:D)
    {
      sigma[d] ~ cauchy(0, ps_sigma);
      mu[d]  ~ normal(0,ps_mu);
    }

  
  for(d in 1:D)
    {
      for(n in 1:N)
        X[n,d] ~ normal(mu_alpha+mu[d], sigma[d]);
    }

}
