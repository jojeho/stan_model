
data{
  int AT;
  int N;
  array[AT] real y;
  array[N] int GT;  
}

parameters{
  real beta;
  real<upper=0> h_alpha;
  real<lower=0> alpha_tau;
  vector[N] alpha_z;
  real alpha1;
  vector<lower=0> sigma;
}


transformed parameters{
  array[N] real alpha;
  alpha[1]=alpha1;
  for( n in 2:N)
    {
      alpha[n]=alpha[n-1]*h_alpha +alpha_tau*alpha_z[n];
    }
}

model{

  alpha_z ~ std_normal();
  alpha_tau ~ cauchy(0,2);
  alpha1 ~ normal(1,5);
  h_alpha~normal(0,5);

  beta ~normal(0,5);
  sigma ~ normal(0,3);

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

generated quantities{
  //
  array[AT-N] real y_hat=rep_array(0,AT-N);
  array[AT-N] real log_lik;
  {
    int pos=1;
    for(n in 1:N)
      {
        int T=GT[n];
        array[T] real yy=segment(y,pos,T);
        for(t in 2:T)
          {
            y_hat[pos+t-1-n] = normal_rng(yy[t-1]*beta+alpha[n] , sigma[n]);
            log_lik[pos+t-1-n]=normal_lpdf(yy[t]|yy[t-1]*beta+alpha[n] , sigma[n]);
          }
        pos +=T;
      }
  }
    
}
  

