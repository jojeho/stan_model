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
  vector<lower=0>[N] sigma;
  vector<lower=0>[N] g_sigma;
}

transformed parameters{
  
  array[N] matrix[T,T] gp_obj ;
  array[N] real zmu;
  array[N] vector[T] mu;
  for(n in 1:N)
    gp_obj[n]= gp(x,tau,smooth,sigma[n]);

  for(n in 1:N)
    {
      mu[n] = gp_obj[n] * eta[n];
    }

}

model{

  tau ~ cauchy(0,Ps_tau);
  smooth~normal(0,Ps_smooth);
  for(n in 1:N)
    eta[n] ~ normal(0,1);
  
  sigma ~normal(0,1);
  g_sigma ~normal(0,2);

  for(n in 1:N)
    {
      for(t in 1:T)
        y[n][t,] ~normal(mu[n][t] ,g_sigma[n]);
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



