functions {
  #include "exp_lib.stan"
}

data{
  int N;
  int T;
  int unit;
  array[N] matrix[T,unit]  y;
  real Ps_tau;
  real Ps_smooth;
}

transformed data{
  
  vector[T] zero_mu;
  array[T] real x;
  for(i in 1:T)
    {
      x[i]=(i-1);
    }
  
  real delta=1e-10;
  zero_mu =rep_vector(0,T);
}


parameters{
  
  real<lower=0> tau;
  real<lower=0> smooth;
  array[N] vector[T] eta;
  real<lower=0> sigma1;

  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  vector<lower=0>[N] sigma;
}

transformed parameters{
  
  matrix[T,T] gp_obj ;
  array[N] vector[T] g_sigma;
  array[N] real zmu;
  array[N] vector[T] mu;
  gp_obj= gp(x,tau,smooth,sigma[n]);

  for(n in 1:N)
    {
      mu[n] = gp_obj * eta[n];
    }

  for(n in 1:N)
    {
      zmu[n]=mean(mu[n]);
    }
  
  for(n in 1:N)
    {
      g_sigma[n][1]=sigma1;
      for (t in 2:T)
        g_sigma[n][t] = sqrt(alpha0
                             + alpha1 * pow(mu[n][t-1] - zmu[n], 2)
                             + beta1 * pow(g_sigma[n][t-1], 2));  
    }
  
}

model{

  tau ~ cauchy(0,Ps_tau);
  smooth~normal(0,Ps_smooth);
  for(n in 1:N)
    eta[n] ~ normal(0,1);
  
  sigma1~ normal(0,1);
  sigma ~normal(0,1);

  for(n in 1:N)
    {
      for(t in 1:T)
        y[n][t,] ~normal(mu[n][t] ,g_sigma[n][t]);
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


