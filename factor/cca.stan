//Canonical Correlation Analysis and Group Factor Analysis
//file:///C:/doc/Probabilistic%20Factor%20Analysis%20Methods.html
data {
  int<lower=0> N; // Number of samples
  int<lower=0> D; // Sum of the original dimensions of every view
  int<lower=0> M; // The number of views
  int<lower=0> K; // The latent dimension
  int<lower = 0> Dm[M]; //The original dimension for each view

  matrix[N, D] X; // The data matrix
}

parameters {
  matrix[N, K] Z; // The latent matrix
  matrix[K, D] W; // The weight matrix. You will be learning *one* single W with structured sparsity - both view-specific components and shared components are captured by this
  vector<lower=0>[M] tau; // View-specific noise terms
  matrix<lower=0>[M,K] alpha; // View-specific ARD prior
}

transformed parameters{
  vector<lower=0>[M] t_tau;
  matrix<lower=0>[M,K] t_alpha; 
  t_alpha = inv(sqrt(alpha));
  t_tau = inv(sqrt(tau));
}

model {
  // You will need to loop through 1:Dm[m] for each view m. Increment ind seperately to index the concatenated X and W.
  int ind;
  tau ~ gamma(1,1);			
  to_vector(Z) ~ normal(0,1); // because sampling K dimensional standard normal multivariate is equivalent to sampling from k univariate standard normal distributions.
  to_vector(alpha) ~ gamma(1e-3,1e-3); // stack columns of alpha to a vector and sample quickly.	
  ind = 0;
  // There is a more efficient way to do this with ragged arrays, but the effort spent is hugely disproportionate to the speed-ups obtained.
  for (m in 1:M) {	
    for (d in 1:Dm[m]) {
      ind = ind + 1;      
      W[,ind] ~ normal(0.0, t_alpha[m,]);
      X[,ind] ~ normal(Z*W[,ind], t_tau[m]);  
    }
  }
