data{
  int AT;
  int N;
  array[AT] real y;
  array[N] int GT;  
}

transformed  data{
  vector[N] h_att;
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      array[T] real yy=segment(y,pos,T);
      h_att[n]=mean(yy);
      pos +=T;
    }

}

parameters{
  array[N] real alpha;
  array[N] real<lower=0> sigma;
  real beta1;
}

transformed parameters{
  array[N] real beta;
  for(n in 2:N)
    {
      beta[n]=beta1*beta[n-1];
    }
}
model{
  for(n in 1:N)
    {
      beta1 ~ normal(0,5);
      alpha[n] ~normal(0,0.1);
    
    }

  sigma ~normal(0,0.1);
  
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      array[T] real yy=segment(y,pos,T);
      for(t in 2:T)
        {
          yy[t] ~ normal(yy[t-1]*beta[n]+alpha[n] , sigma[n]);
        }
      pos +=T;
    }
}

/* generated quantities{ */
/*   array[AT] real y_hat;//=rep_array(0,AT); */
/*   array[AT] real log_lik; */
/*   int pos=1; */
/*   for(n in 1:N) */
/*     { */
/*       int T=GT[n]; */
/*       array[T] real yy=segment(y,pos,T); */
/*       for(t in 2:T) */
/*         y_hat[pos+t-1] = normal_rng(yy[t-1]*beta[n]+alpha[n] , sigma[n]); */
      
/*       pos +=T; */
/*     } */
/* } */
  

