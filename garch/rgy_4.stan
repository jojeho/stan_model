
data{
  #include "../input/rgy.stan"
}

transformed data {
  real alpha_prior;
  real sigma_prior;
  real mu_prior;
  mu_prior=0.1;
  alpha_prior=0.1;
  sigma_prior=1;
    
}


parameters {

  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  //array[N] real<lower=0> sigma1;


  //real mu_tau;
  //vector[N] mu_z;
  //real mu_raw;

  //real sigma1_tau;
  //vector[N] sigma1_z;
  //real sigma1_raw;
  real sigma1;
  vector[N] mu;
}

transformed  parameters{
  //array[N] real  mu;
  /* for(n in 1:N) */
  /*   { */
  /*     mu[n]= mu_raw+ mu_tau*mu_z[n]; */
  /*   } */
}

model {

  alpha1 ~ normal(0,1);
  alpha0 ~ normal(0,alpha_prior);
  beta1 ~ normal(0,1);
  sigma1 ~ normal(0,sigma_prior);

  //mu_z ~ std_normal();
  //mu_tau ~ cauchy(0,3);
  mu ~ normal(0,0.1);


  sigma1 ~ normal(0,0.1);  
  
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      vector[T] yy=segment(y,pos,T);
      array[T] real sigma;
      sigma[1]=sigma1;
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

