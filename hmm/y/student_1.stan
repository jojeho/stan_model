data{
  int T;
  int K;
  real y[T];
  real dir_tr;
  real dir_rho;
  int v;
  real mu_scale_prior;
}

parameters{
  ordered [K] mu;
  //vector<lower=0>[K] sigma;
  //real<lower=0> sigma;
  
  simplex[K] tr[K];
  simplex[K] rho;

  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  real<lower=0> sigma1;
}

transformed parameters{
  matrix[K,K] Gamma= rep_matrix(0, K, K);
  vector<lower=0>[T] sigma;
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }

  sigma[1] = sigma1;
  for (t in 2:T)
    sigma[t] = sqrt(
                     alpha1 * pow(y[t-1] , 2)
                    + beta1 * pow(sigma[t-1], 2));
  
}


model{

  sigma1 ~ normal(0,0.1);
  alpha1 ~ normal(0,1);
  beta1 ~ normal(0,5);
  //sigma ~ cauchy(0,5);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho  ~dirichlet(rep_vector(dir_rho/K,K));
  mu ~normal(0,mu_scale_prior);

  
  matrix[K,T] ob;

  
  for(k in 1:K)
    {
      for(t in 1:T)
        {
          //ob[k,t]=student_t_lpdf(y[t]|v,mu[k] ,sigma);
          ob[k,t]=normal_lpdf(y[t]|mu[k] ,sigma[t]);
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
          ob[k,t]=normal_lpdf(y[t]|mu[k] ,sigma[t]);
        }
    }
  
  matrix[K,T] prob = hmm_hidden_state_prob(ob,Gamma,rho);
  
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
