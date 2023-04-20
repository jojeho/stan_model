data{
  #include "../input/rgy.stan"
}


parameters {
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  //array[N] real<lower=0> sigma1;
  real<lower=0> sigma1;
  
}

model {

  sigma1 ~normal(0,2);
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
                          + alpha1 * pow(yy[t-1] , 2)
                          + beta1 * pow(sigma[t-1], 2));
        }
      
      yy~ normal(0, sigma);      
      
      pos +=T;
    }
  
}


