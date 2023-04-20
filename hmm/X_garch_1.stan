
data{
  int N;
  int D;
  array[N] vector[D] X;
  int K;
  real dir_tr;
  real dir_rho;
  //  real<lower=2> v;
}


parameters{
  //vector<lower=0>[K] sigma;
  //  real<lower=0> sigma;
  simplex[K] tr[K];
  simplex[K] rho;

  ordered[K] mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=0> sigma1;

  real<lower=0> s_beta;

}

transformed parameters{
  real<lower=0> sigma[N];
  matrix[K,K] Gamma= rep_matrix(0, K, K);
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }
  matrix[K,D] ob=rep_matrix(0,K,D);
  for(d in 1:D)
    {
      for(k in 1:K)
        {

          sigma[1] = sigma1;
          for (n in 2:N)
            sigma[n] = sqrt(alpha0
                            + alpha1 * pow(X[n-1,d] - mu[k], 2)
                            + beta1 * pow(sigma[n-1], 2));          
          for(n in 2:N)
            {
              ob[k,d] +=normal_lpdf(X[n,d]|mu[k],sigma);
            }
        }
    }
  
}


model{
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  mu ~ normal(0,0.5);
  alpha0~normal(0,1);
  sigma1 ~normal(0,3);
  s_beta ~normal(0,1);
  
  
  target += hmm_marginal(ob,Gamma,rho);
}

generated quantities{
  matrix[K,D] prob=hmm_hidden_state_prob(ob,Gamma,rho);
}
