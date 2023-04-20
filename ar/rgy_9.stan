
data{
  int AT;
  int N;
  array[AT] real y;
  array[N] int GT;  
}

parameters{
  array[N] real alpha;
  real beta;
  vector<lower=0>[N] sigma;
}

model{
  beta ~ normal(0,5);
  for(n in 1:N)
    {
      alpha[n] ~normal(0,1);
    }

  sigma ~normal(0,4);
  
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      array[T] real yy=segment(y,pos,T);
      for(t in 2:T)
        {
          yy[t] ~ normal(yy[t-1]*beta+alpha[n] , sigma[n]);
        }
      pos +=T;
    }
}

/* generated quantities{ */
/*   // */
/*   array[AT-N] real y_hat=rep_array(0,AT-N); */
/*   array[AT-N] real log_lik; */
/*   { */
/*     int pos=1; */
/*     for(n in 1:N) */
/*       { */
/*         int T=GT[n]; */
/*         array[T] real yy=segment(y,pos,T); */
/*         for(t in 2:T) */
/*           { */
/*             y_hat[pos+t-1-n] = normal_rng(yy[t-1]*beta+alpha[n] , sigma[n]); */
/*             log_lik[pos+t-1-n]=normal_lpdf(yy[t]|yy[t-1]*beta+alpha[n] , sigma[n]); */
/*           } */
/*         pos +=T; */
/*       } */
/*   } */
/* } */
  



