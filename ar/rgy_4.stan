
data{
  int AT;
  int N;
  array[AT] real y;
  array[N] int GT;  
}

parameters{
  real beta;
  real alpha1;
  real<lower=0> sigma1;
  real<upper=0> h_alpha;
  real<lower=0> h_sigma;

  real<lower=0> alpha_tau;
  vector[N] alpha_z;

  real<lower=0> sigma_tau;
  vector[N] sigma_z;
}


transformed parameters{
  array[N] real alpha;
  array[N] real<lower=0> sigma;
  alpha[1]=alpha1;
  sigma[1]=sigma1;
  for( n in 2:N)
    {
      alpha[n]=alpha[n-1]*h_alpha +alpha_tau*alpha_z[n];
      sigma[n]=sigma[n-1]*h_sigma +sigma_tau*sigma_z[n];
    }
}

model{

  alpha_z ~ std_normal();
  alpha_tau ~ cauchy(0,2);

  sigma_tau ~ cauchy(0,2);
  sigma_z ~ std_normal();
  
  alpha1 ~ normal(1,5);
  beta ~normal(0,5);
  sigma1~cauchy(0,4);
  h_alpha~normal(0,5);
  h_sigma~normal(0,5);

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
  

