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
  int unit;
  int K;
  real dir_tr;
  real dir_rho;
  
  array[T] vector[unit]  y;
}


parameters{
  
  real<lower=0.001> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=0> sigma1;

  ordered[K] mu;
  simplex[K] tr[K];
  simplex[K] rho;
}

transformed parameters{

  matrix[K,K] Gamma= rep_matrix(0, K, K);
  real mm = mean(mu);
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }

    
  
  array[T]  vector<lower=0>[unit] sigmas;
  for(t in 1:T)
    {
      sigmas[t,1] = sigma1;
      for (n in 2:unit)
        sigmas[t,n] = sqrt(alpha0
                           + alpha1 * pow(y[t,n-1] - mm, 2)
                           + beta1 * pow(sigmas[t,n-1], 2));

    }
  
}



model{

  //tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  //rho ~dirichlet(rep_vector(dir_rho/K,K));
  //vector[2] v=[6,6]';
  //tr  ~ dirichlet(v);
  mu  ~normal(0,1);
  alpha0~normal(0,1);
  sigma1 ~normal(0,1);
  
  matrix[K,T] ob;      
  for(k in 1:K)
    {
      for(t in 1:T)
        {
          ob[k,t]=normal_lpdf(y[t,]|mu[k] ,sigmas[t,]);
        }
    }
  
  target += hmm_marginal(ob,Gamma,rho);

  //  #include "normal_loglik.lstan"
}

generated quantities{
  matrix[K,T] prob;
  matrix[K,T] ob;

  for(k in 1:K)
    {
      for(t in 1:T)
        {
          ob[k,t]=normal_lpdf(y[t,]|mu[k] ,sigmas[t,]);
        }
    }

  prob=hmm_forward_prob(ob,Gamma,rho); 
}




//#include "normal_gen.lstan"

