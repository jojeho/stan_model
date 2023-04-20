
data{
  int AT;
  int N;
  array[AT] real y;
  array[N] int GT;

  int NT;
  array[NT] real y_n;
  
}

parameters{
  array[N] real alpha;
  array[N] real beta;
  array[N] real<lower=0> sigma;
}

/* model{ */
/*   for(n in 1:N) */
/*     { */
/*       beta[n] ~ normal(0,1); */
/*       alpha[n] ~normal(0,1); */
/*     } */

/*   sigma ~normal(0,1); */
  
/*   int pos=1; */
/*   for(n in 1:N) */
/*     { */
/*       int T=GT[n]; */
/*       array[T] real yy=segment(y,pos,T); */
/*       for(t in 2:T) */
/*         { */
/*           yy[t] ~ normal(yy[t-1]*beta[n]+alpha[n] , sigma[n]); */
/*         } */
/*       pos +=T; */
/*     } */
/* } */

generated quantities{
  //
  array[AT-N] real y_hat=rep_array(0,AT-N);
  array[AT-N] real log_lik;
  array[N] vector[NT-1] y_hat_n;
  {
    int pos=1;
    for(n in 1:N)
      {
        int T=GT[n];
        array[T] real yy=segment(y,pos,T);
        for(t in 2:T)
          {
            //y_hat[pos+t-1-n] = normal_rng(yy[t-1]*beta[n]+alpha[n] , sigma[n]);
            log_lik[pos+t-1-n]=normal_lpdf(yy[t]|yy[t-1]*beta[n]+alpha[n] , sigma[n]);
          }
        pos +=T;
      }

    for(n in 1:N)
      {
        for(t in 2:NT)
          {
            y_hat_n[n,t-1]=normal_rng(y_n[t-1]*beta[n]+alpha[n],sigma[n]);
          }
      }
  }
}
  

