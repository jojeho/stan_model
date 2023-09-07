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
  array[T] vector[unit]  y;
}


parameters{
  
  real<lower=0.001> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=0> sigma1;
  real mu1;
  real<lower=0> mu_alpha;
  real mu_beta;
  vector[T] eta;
}

transformed parameters{

  array[T]  vector[unit] sigmas;
  vector[T]  mu ;
  mu[1]=mu1;
  
  for(t in 2:T)
    {
      mu[t]= mu[t-1]*mu_alpha+mu_beta*eta[t];
    }
  
  for(t in 1:T)
    {
      sigmas[t,1] = sigma1;
      for (n in 2:unit)
        sigmas[t,n] = sqrt(alpha0
                           + alpha1 * pow(y[t,n-1] - mu[t], 2)
                           + beta1 * pow(sigmas[t,n-1], 2));

    }
  
}



model{

  eta  ~ std_normal();
  mu_beta  ~ normal(0,3);
  alpha0~normal(0,1);
  mu1 ~normal(0,0.1);
  mu_alpha ~ normal(0,1);
  sigma1 ~normal(0,3);  

  for(t in 1:T)
    {
      y[t,] ~ normal(mu[t],sigmas[t,]);
    }
  

  //  #include "normal_loglik.lstan"
}

//#include "normal_gen.lstan"


/* generated quantities{ */
/*   matrix[K,AT] prob; */
/*   { */
/*   int pos=1; */
/*   for(n in 1:N) */
/*     { */
/*       int T=GT[n]; */
/*       array[T] real yy=segment(y,pos,T); */
/*       int nsize=T/unit; */
/*       matrix[K,nsize] ob; */
/*       matrix[K,nsize] local_prob=rep_matrix(0,K,nsize); */
/*       for(k in 1:K) */
/*         { */
/*           int pos2=1; */
/*           for(ni in 1:nsize) */
/*             { */
/*               array[unit] real y_n=segment(y,pos+pos2,unit); */
/*               ob[k,ni]=normal_lpdf(y_n|mu ,sigma[k]); */
/*             } */
/*           pos2=pos2+unit; */
/*         } */


/*       local_prob = hmm_forward_prob(ob,Gamma,rho); */
/*        for(ni in 1:nsize) */
/*           { */
/*             prob[:,pos+ni-1]= local_prob[:,ni]; */
/*           } */
      
/*       pos +=T; */
      
/*     } */
/*   } */
/* } */


