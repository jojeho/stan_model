
data{
  int AT;
  int N;
  array[AT] real y;
  array[N] int GT;
}

parameters{
  real alpha;
  real beta;
  //array[N] real<lower=0> sigma;
  real<lower=0> sigma;
}

model{
  beta ~ normal(0,5);
  alpha ~normal(0,1);

  sigma ~normal(0,3);
  
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      array[T] real yy=segment(y,pos,T);
      for(t in 2:T)
        {
          yy[t] ~ normal(yy[t-1]*beta+alpha , sigma);
        }
      pos +=T;
    }
}

generated quantities{
  array[AT] real y_hat=rep_array(0,AT);
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      array[T] real yy=segment(y,pos,T);
      for(t in 2:T)
        y_hat[pos+t-1] = normal_rng(yy[t-1]*beta+alpha , sigma);
      pos +=T;
    }
}
  

