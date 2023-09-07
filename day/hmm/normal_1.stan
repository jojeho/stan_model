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
  int unit;
}


parameters{

  positive_ordered[K] sigma;
  simplex[K] tr[K];
  simplex[K] rho;
  real mu;
}

transformed parameters{

  matrix[K,K] Gamma= rep_matrix(0, K, K);
  for(k in 1:K)
    {
      Gamma[k,] = tr[k]';
    }
}


model{
  sigma ~ normal(0,2);
  tr  ~ dirichlet(rep_vector(dir_tr/(K),K));
  rho ~dirichlet(rep_vector(dir_rho/K,K));
  mu  ~normal(0,1);

  int pos=1;
  for(n in 1:N)
    {
      int AN=GT[n];
      array[AN] real yy=segment(y,pos,AN);
      int T=AN/unit;
      matrix[K,T] ob;      
      for(k in 1:K)
        {
          int pos2=1;
          for(t in 1:T)
            {
              array[unit] real y_n=segment(yy,pos2,unit);
              ob[k,t]=normal_lpdf(y_n|mu ,sigma[k]);
              pos2=pos2+unit;
            }
        }
      target += hmm_marginal(ob,Gamma,rho);
      pos =pos + AN;
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


