data{
  #include "../../input/X.stan"
}


transformed data {
  real log_unif;
  log_unif = -log(N);
}
parameters {
  real<lower=0> sigma;

  real mu_tilde1;
  real mu_tilde2;
  //  real mu_c;
  real tau;
}

transformed parameters {
  matrix[N,D] lp;
  real mu1;
  real mu2;

  mu1 =  mu_tilde1*tau;
  mu2 =  mu_tilde2*tau;
  
  for(d in 1:D)
    {
      lp[,d] = rep_vector(log_unif, N);
      for (s in 1:N) {
        for (t in 1:N) {
          lp[s,d] = lp[s,d] + normal_lpdf(X[t,d] | t < s ? mu1 : mu2 , sigma );
        }
      }
    }
}

model {

  mu_tilde1 ~ normal(0,1);
  mu_tilde2 ~ normal(0,1);
  tau ~ cauchy(0,5);
  //  mu_c ~ normal(0,1);

  //mu1 ~ normal(0,ps_mu);
  //mu2 ~ normal(0,ps_mu);
  
  sigma ~ normal(0,5);
  //  sigma2 ~ normal(0,ps_sigma);

  for(d in 1:D)
    {
      target += log_sum_exp(lp[,d]);
    }
}

generated quantities {
  vector<lower=1, upper=N>[D] s;
  for(d in 1:D)
    s[d] = categorical_logit_rng(lp[,d]);
}

