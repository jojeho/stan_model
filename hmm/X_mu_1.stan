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
  int D;
  array[N] vector[D] X;
  real dir_tr;
  real dir_rho;
  //  real sigma_prior;    
  
}

transformed data{
  int K =2;
  vector[N] x;
  for(t in 1:N)
    {
      x[t]=t*0.1;
    }  
}

parameters{
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  simplex[K] tr[K];
  //  simplex[K] rho;
  real<lower=0> mu1;
  real<upper=0> mu2;
  //  real<lower=0> kappa;
  
}

transformed parameters{
  vector[K] mu;
  simplex[K] rho;
  mu[1]=mu1;
  mu[2]=mu2;
  rho[1]=0.5;
  rho[2]=0.5;
  
  
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
          for(n in  1:N)
            {
              ob[k,d] +=normal_lpdf(X[n,d]|mu[k] ,sigma);
            }
          ob[k,d] -=log(N);
        }
    }
  
}


model{
  sigma ~ normal(0,2);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  /* beta1 ~ von_mises(pi()/4,kappa); */
  /* beta2 ~ von_mises(-pi()/4,kappa); */
  mu1 ~ normal(-3,1);
  mu2 ~ normal(5.5,1);  
  //kappa~normal(5,2);
  target += hmm_marginal(ob,Gamma,rho);
}

generated quantities{
  //matrix[K,D] prob2=hmm_hidden_state_prob(ob,Gamma,rho);
  matrix[K,D] prob=hmm_forward_prob(ob,Gamma,rho);
}
