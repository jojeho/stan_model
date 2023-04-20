data{
  #include "../input/y.stan"
}


parameters {
  real alpha0;
  real alpha1;
  real beta1;
  //array[N] real<lower=0> sigma1;
  real<lower=0> sigma1;
  
}

generated quantities{
  array[T] real sigma;
  sigma[1]=sigma1;
  for (t in 2:T)
    {
      sigma[t] = sqrt(alpha0
                      + alpha1 * pow(y[t-1] , 2)
                      + beta1 * pow(sigma[t-1], 2));
    }
}



