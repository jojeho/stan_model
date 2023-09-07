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
  int T;
  int K;
  real y[T];
  real dir_tr;
  real dir_rho;
}

parameters{
  positive_ordered [K] mu;
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  
  simplex[K] tr[K];
  simplex[K] rho;
  real<lower=0> v;
}

transformed parameters{
  matrix[K,T] ob;
    
  for(k in 1:K)
    {
      /* for(t in 1:T) */
      /*   { */
      /*     ob[k,t]=normal_lpdf(y[t]|mu[k] ,sigma); */
      /*   } */
      for(t in 1:T)
      {
        ob[1,t]=student_t_lpdf(y[t]|v,mu[1],sigma);
        ob[2,t]=student_t_lpdf(y[t]|v,mu[2],sigma);
      }
        
    }

  matrix[K,K] Gamma= rep_matrix(0, K, K);

  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }

  }


model{
  
  //sigma[1] ~ normal(0,1);
  sigma ~ normal(0,1);
  //sigma[2] ~ normal(0,1);
  v ~ gamma(1,0.1);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  mu[1] ~ normal(0.2,1);
  mu[2] ~ normal(0.5,2);
  //mu ~normal(0,1);

  target += hmm_marginal(ob,Gamma,rho);
}

generated quantities{
  matrix[K,T] prob=hmm_forward_prob(ob,Gamma,rho);
 /*  array[T] int label=hmm_latent_rng(ob,Gamma,rho); */
 /*  vector[T] y_hat; */
 /* stud vector[T] log_lik; */
 /*  for(t in 1:T) */
 /*    { */
 /*      y_hat[t]=normal_rng(mu[label[t]],sigma); */
 /*      log_lik[t]=normal_lpdf(y[t]|mu[label[t]],sigma); */
 /*    } */
}



