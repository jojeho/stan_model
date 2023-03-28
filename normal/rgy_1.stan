data{
  #include "../input/rgy.stan"
}

transformed data{
}

parameters {
  real<lower=0> mu_tau;
  real mu_z;
  real mu_beta;

  real<lower=0> sigma_tau;
  real sigma_z;
  real sigma_beta;

  real mu_init;
  real<lower=0> sigma_init;
}

transformed parameters{
  vector[N] mu;
  vector[N] sigma;

  mu[1]=mu_init;
  sigma[1]=sigma_init;
    
  for(n in 2:N)
    {
      mu[n]= mu[n-1]*mu_beta+ mu_tau*mu_z[n];
      sigma[n]=sigma[n-1]*sigma_beta + sigma_tau*sigma_z[n];
    }
}

model {

  mu_init~normal(0,3);
  sigma_init~normal(0,3);
  
  sigma_z ~ std_normal();
  sigma_tau ~ cauchy(0,3);
  sigma_beta ~ normal(0,1);

  mu_z ~ std_normal();
  mu_tau ~ cauchy(0,3);
  mu_beta ~ normal(0,1);  

  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      vector[T] yy=segment(y,pos,T);
      array[T] real sigma;
      yy~ normal(mu[n], sigma[n]);      
      pos +=T;
    }
  
}

