data{
  #include "../input/rgy.stan"
}

transformed data{
}

parameters {
  array[N] real  mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  //array[N] real<lower=0> sigma1;
  vector<lower=0>[N] sigma1;
  
}

model {

  sigma1 ~normal(0,5);
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

