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
  real v;
  real mu_scale_prior;
}

parameters{
  ordered [K] mu;
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

}


model{

  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  mu ~normal(0,mu_scale_prior);
  sigma ~ normal(0,2);
  matrix[K,T] ob;

  
  for(k in 1:K)
    {
      for(t in 1:T)
        {
          ob[k,t]=student_t_lpdf(y[t]|v,mu[k] ,sigma);
        }
    }

  
  target += hmm_marginal(ob,Gamma,rho);
}

generated quantities{
  matrix[K,T] ob;
  for(k in 1:K)
    {
      for(t in 1:T)
        {
          ob[k,t]=student_t_lpdf(y[t]|v,mu[k] ,sigma);
        }
    }
  
  matrix[K,T] prob = hmm_forward_prob(ob,Gamma,rho);
  
}
/*   vector[K] n_prob=rep_vector(0,K); */
/*   vector[T] log_lik; */

/*   for(k in 1:K) */
/*     { */
/*       for(j in 1:K) */
/*         { */
/*           n_prob[k] +=log(prob[j,T])+log(Gamma[j,k]); */
/*         } */
/*     } */
/*   n_prob=softmax(n_prob); */
  
/*   for(t in 1:T) */
/*     { */
/*       log_lik[t]=log_sum_exp(log(prob[,t])+ob[,t]); */
/*     } */
/* } */
