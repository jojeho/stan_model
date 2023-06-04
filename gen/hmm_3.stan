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


#include "input.lstan"

transformed data{
  int K=4;
}

parameters{
  //vector<lower=0>[K] sigma;
  real<lower=0> sigma;
  
  simplex[K] tr[K];
  simplex[K] rho;

  real mu1;
  real mu2;
  real mu3;
  real mu4;
}

transformed parameters{
  matrix[K,T] ob;
  vector[K] mu;
  mu[1]=mu1;
  mu[2]=mu2;
  mu[3]=mu3;
  mu[4]=mu4;
  
  for(k in 1:K)
    {
      for(t in 1:T)
        {
          ob[k,t]=normal_lpdf(y[t]|mu[k] ,sigma);
        }
    }

  matrix[K,K] Gamma= rep_matrix(0, K, K);

  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }


  }


model{
  sigma ~ normal(0,3);
  tr  ~ dirichlet(rep_vector(1,K));
  rho  ~dirichlet(rep_vector(1,K));
  mu1 ~normal(-2,2);
  mu2 ~normal(-1,2);
  mu3 ~normal(1,2);
  mu4 ~normal(2,2);

  target += hmm_marginal(ob,Gamma,rho);
}

generated quantities{
  vector[NT] y_hat;
  vector[NT] log_lik;
  {
    matrix[K,NT] n_ob;
    
    for(k in 1:K)
      {
        for(t in 1:NT)
          {
            n_ob[k,t]=normal_lpdf(n_y[t]|mu[k] ,sigma);
          }
      }

    matrix[K,NT] prob=hmm_forward_prob(n_ob,Gamma,rho);
    
    real p_mu;
    for(t in 1:NT)
      {
        p_mu = dot_product(mu,prob[,t]);
        y_hat[t]=normal_rng(p_mu,sigma);
        log_lik[t]=normal_lpdf(n_y[t]|p_mu,sigma);
      }
  }
}




