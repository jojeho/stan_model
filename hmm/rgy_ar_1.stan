functions {
  
  // ref: https://github.com/luisdamiano/gsoc17-hhmm/blob/master/hmm/stan/hmm-multinom.stan
  matrix hmm_forward_prob(matrix log_omega, matrix Gamma, vector rho) {
    
    int n_state = rows(log_omega);
    int n_obs = cols(log_omega);
    matrix[n_state, n_obs] unnorm_log_alpha; 
    matrix[n_state, n_obs] alpha; // p(x_t=j|y_{1:t})
    array[n_state] real accumulator;
    
    // first alpha
    unnorm_log_alpha[,1] = log(rho) + log_omega[,1];
    
    // other alphas
    for (t in 2:n_obs) {
      for (j in 1:n_state) {
        for (i in 1:n_state) {
          accumulator[i] = unnorm_log_alpha[i, t-1] + log(Gamma[i, j]) + log_omega[j, t];
        }
        unnorm_log_alpha[j, t] = log_sum_exp(accumulator);
      }
    }

    // normalize alpha later for numerical stability
    for (t in 1:n_obs) {
      alpha[,t] = softmax(unnorm_log_alpha[,t]);
    }
    return alpha;
  } // !hmm_forward_prob
}



data{
  int N;
  int AT;
  array[AT] real  y;
  array[N] int GT;
  int K;
  real dir_tr;
  real dir_rho;
}

parameters{
  ordered[K] alpha;
  real beta;
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  simplex[K] tr[K];
  simplex[K] rho;

}

transformed parameters{
  matrix[K,K] Gamma ;//= rep_matrix(0, K, K);
  for(k in 1:K)
    Gamma[k,]=tr[k]';
  
}

model{
  sigma ~ cauchy(0,3);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  alpha ~normal(0,2);
  beta ~ normal(0,1);

  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];

      matrix[K,T-1] ob;

      array[T] real yy=segment(y,pos,T);
      
  
      for(k in 1:K)
        {
          for(t in 2:T)
            {
              ob[k,t-1]=normal_lpdf(yy[t]|beta*yy[t-1]+alpha[k] ,sigma);
            }
        }

      target += hmm_marginal(ob,Gamma,rho);
      pos +=T;
    }

}



