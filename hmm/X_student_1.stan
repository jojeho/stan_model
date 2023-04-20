
data{
  int N;
  int D;
  array[N] vector[D] X;
  int K;
  real dir_tr;
  real dir_rho;
  real<lower=2> v;
}


parameters{
  ordered[K] mu;
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  simplex[K] tr[K];
  simplex[K] rho;

}

transformed parameters{
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
          for(n in 1:N)
            {
              ob[k,d] +=student_t_lpdf(X[n,d]|v,mu[k] ,sigma);
            }
        }
    }
  
}


model{
  sigma ~ normal(0,1);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  mu ~normal(0,0.5);
  target += hmm_marginal(ob,Gamma,rho);
}

generated quantities{
  matrix[K,D] prob=hmm_hidden_state_prob(ob,Gamma,rho);
}
