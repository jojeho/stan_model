data{
  #include "../input/rgy.stan"
}

transformed data{
}

parameters {
  //  vector[N]  mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  //array[N] real<lower=0> sigma1;

  real<lower=0> mu_tau;
  vector[N] mu_z;
  real mu_beta;

  real<lower=0> sigma1_tau;
  vector[N] sigma1_z;
  real sigma1_beta;

  real mu_init;
  real<lower=0> sigma1_init;
}

transformed parameters{
  vector[N] mu;
  vector[N] sigma1;

  mu[1]=mu_init;
  sigma1[1]=sigma1_init;
    
  for(n in 2:N)
    {
      mu[n]= mu[n-1]*mu_beta+ mu_tau*mu_z[n];
      sigma1[n]=sigma1[n-1]*sigma1_beta + sigma1_tau*sigma1_z[n];
    }
}

model {

  mu_init~normal(0,1);
  sigma1_init~normal(0,1);
  
  sigma1_z ~ std_normal();
  sigma1_tau ~ cauchy(0,3);
  sigma1_beta ~ normal(0,1);

  mu_z ~ std_normal();
  mu_tau ~ cauchy(0,3);
  mu_beta ~ normal(0,1);  

  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      vector[T] yy=segment(y,pos,T);
      array[T] real sigma;
      sigma[1]=sigma1[n];
      for (t in 2:T)
        {
          sigma[t] = sqrt(alpha0
                          + alpha1 * pow(yy[t-1] - mu[n], 2)
                          + beta1 * pow(sigma[t-1], 2));
        }
      
      yy~ normal(mu[n], sigma);      
      
      pos +=T;
    }
  
}

