data{
#include "../input/rgy.stan"
#include "../input/gen_ny_i.stan"
}

transformed data{
  int ANT=NT+1;
}

parameters {
  array[N] real alpha;
  array[N] real beta;
  array[N] real<lower=0> sigma;
}

model{
  beta ~ normal(1,5);
  alpha ~normal(0,0.1);

  sigma ~normal(0,0.1);
  
  int pos=1;
  for(n in 1:N)
    {
      int T=GT[n];
      vector[T] yy=segment(y,pos,T);
      for(t in 2:T)
        {
          yy[t] ~ normal(yy[t-1]*beta[n]+alpha[n] , sigma);
        }
      pos +=T;
    }
}


generated quantities{

  array[ANT] real n_yhat;
  array[N] real last_y;

  int pos=1;

  for(n in 1:N)
    {
      int T=GT[n];
      vector[T]  yy=segment(y,pos,T);
      last_y[n]=yy[T];
      pos +=T;
    }

  
  
  for(t in 1:ANT)
    {

      if( t == 1)
        {
          n_yhat[1]=normal_rng(last_y[index]*beta[index]+alpha[index] , sigma[index]);
          
        }
      else
        {
          n_yhat[t]=normal_rng(n_y[t-1]*beta[index]+alpha[index] , sigma[index]);
        }
    }
  
}



