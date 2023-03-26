data{
 #include "../../input/X.stan"
}

parameters{
  
  real mu_alpha;
  vector<lower=0>[D] sigma;
  //real<lower=0> ps_mu;
  real<lower=0> ps_sigma;
  real tau;
  vector[D] tild;
}

transformed parameters{
  vector[D] mu;
  for(d in 1:D)
    {
      mu[d]=tild[d]*tau +mu_alpha;
    }
}
model{

  #ps_mu~normal(0,1);
  //mu ~normal(0,0.05);
  ps_sigma ~ normal(0,5);
  //ps_mu ~ normal(0,5);
  mu_alpha ~normal(0,0.1);  
  tau ~ cauchy(0,10);
  tild ~ normal(0,1);
  for(d in 1:D)
    {
      
      sigma[d] ~ cauchy(0, ps_sigma);
      //mu[d]  ~ normal(0,ps_mu);
    }

  
  for(d in 1:D)
    {
      for(n in 1:N)
        X[n,d] ~ normal(mu_alpha+mu[d], sigma[d]);
    }

}
