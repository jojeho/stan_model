
data{
  int AT;
  int N;
  array[AT] real y;
  array[N] int GT;  
}

parameters{
  real alpha;
  real beta1;
  real<lower=0> sigma1;
  real h_beta;
  real<lower=0> h_sigma;
}

transformed parameters{
  array[N] real beta;
  array[N] real<lower=0> sigma;
  beta[1]=beta1;
  sigma[1]=sigma1;
  for( n in 2:N)
    {
      beta[n]=beta[n-1]*h_beta;
      sigma[n]=sigma[n-1]*h_sigma;
    }
}

model{

  beta1 ~ normal(1,5);
  alpha ~normal(0,0.1);
  sigma1~cauchy(0,4);
  h_beta~normal(0,5);
  h_sigma~normal(0,5);

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

generated quantities{
  //
  array[AT-N] real y_hat=rep_array(0,AT-N);
  array[AT-N] real log_lik;
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      array[T] real yy=segment(y,pos,T);
      for(t in 2:T)
        {
          y_hat[pos+t-1-n] = normal_rng(yy[t-1]*beta[n]+alpha[n] , sigma[n]);
          log_lik[pos+t-1-n]=normal_lpdf(yy[t]|yy[t-1]*beta[n]+alpha[n] , sigma[n]);
        }
      
      
      pos +=T;
    }
}
  

