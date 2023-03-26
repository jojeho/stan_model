data{
  #include "../../input/X.stan"
}


transformed data {
  real log_unif;
  log_unif = -log(N);
}
parameters {
  real mu;
  vector<lower=0>[D] sigma1;
  vector<lower=0>[D] sigma2;
  //vector<lower=0>[D] sigma2;
  real<lower=0> ps_mu;
  real<lower=0> ps_sigma; 
}

transformed parameters {
  matrix[N,D] lp;
  
  
  for(d in 1:D)
    {
      lp[,d] = rep_vector(log_unif, N);
      for (s in 1:N) {
        for (t in 1:N) {
          lp[s,d] = lp[s,d] + normal_lpdf(X[t,d] | mu , t < s ? sigma1[d]:sigma2[d] );
        }
      }
    }
  
}
model {

  ps_mu ~ normal(0,10);
  ps_sigma ~ normal(0,10);
  mu ~ normal(0,ps_mu);
  sigma1 ~ normal(0,ps_sigma);
  sigma2 ~ normal(0,ps_sigma);
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
