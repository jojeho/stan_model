data{
  #include "../input/rgyi.stan"
}

transformed data{
  int ANT=NT+1;

  /* { */
  /*   int pos=1; */
  /*   for(n in 1:N) */
  /*     { */
  /*       int T=GT[n]; */
  /*       array[T] real yy=segment(y,pos,T); */
  /*       last_y[n]=yy[T]; */
  /*       pos +=T; */
  /*     } */
  /* } */

}

parameters {
  real mu;
  real<lower=0> alpha0;
  real<lower=0,upper=1> alpha1;
  real<lower=0,upper=(1-alpha1)> beta1;
  //array[N] real<lower=0> sigma1;
  real<lower=0> sigma1;
  
}
transformed parameters {
  

}
model {

  sigma1 ~normal(0,0.2);
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
                          + alpha1 * pow(yy[t-1] - mu, 2)
                          + beta1 * pow(sigma[t-1], 2));
        }

      yy~ normal(mu, sigma);      
      
      pos +=T;
    }
  
}

generated quantities{

  array[ANT] real n_sigma;

  {
    array[N] real last_sigma;
    array[N] real last_y;
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
                            + alpha1 * pow(yy[t-1] - mu, 2)
                            + beta1 * pow(sigma[t-1], 2));
          }
        last_sigma[n]=sigma[T];
        last_y[n]=yy[T];
        pos +=T;
      }

  
  
    for(t in 1:ANT)
      {

        if( t == 1)
          {
            n_sigma[1]=sqrt(alpha0
                            + alpha1 * pow(last_y[index] - mu, 2)
                            + beta1 * pow(last_sigma[index], 2));
          
          }
        else
          {
            n_sigma[t]=sqrt(alpha0
                            + alpha1 * pow(n_y[t-1] - mu, 2)
                            + beta1 * pow(n_sigma[t-1], 2));
          }
      }

  }
}

