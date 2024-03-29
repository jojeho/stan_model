//Principal Component Analysis
//file:///C:/doc/Probabilistic%20Factor%20Analysis%20Methods.html
data {
/* int<lower=0> N; // Number of samples */
/* int<lower=0> D; // The original dimension */
/*  matrix[N, D] X; // The data matrix */
  
#include "../input/X.stan"
int<lower=0> K; // The latent dimension

}

parameters {
  matrix[N, K] Z; // The latent matrix
  matrix[D, K] W; // The weight matrix
  real<lower=0> tau; // Noise term 
  vector<lower=0>[K] alpha; // ARD prior
}

transformed parameters{
  vector<lower=0>[K] t_alpha;
  real<lower=0> t_tau;

  t_alpha = inv(sqrt(alpha));
  t_tau = inv(sqrt(tau));
}

model {
  tau ~ gamma(1,1);			
  to_vector(Z) ~ normal(0,1);
  alpha ~ gamma(1e-3,1e-3);
  //alpha ~ gamma(20,20);				
  for(k in 1:K) W[,k] ~ normal(0, t_alpha[k]);
  to_vector(X) ~ normal(to_vector(Z*W'), t_tau);

}
