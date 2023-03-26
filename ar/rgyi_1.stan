data{
  #include "../input/rgyi.stan"
    }

transformed data{
  int ANT=NT+1;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
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
          yy[t] ~ normal(yy[t-1]*beta+alpha , sigma);
        }
      pos +=T;
    }
}


generated quantities{

  array[ANT] real n_yhat;
  array[N] real last_yhat;
  array[N] real last_y;

  int pos=1;

  for(n in 1:N)
    {
      int T=GT[n];
      vector[T] yy=segment(y,pos,T);
      for (t in 2:T)
        {
          n_yhat[t] =normal_rng(yy[t-1]*beta+alpha , sigma);
          
        }
      last_y[n]=yy[T];
      pos +=T;
    }

  
  
  for(t in 1:ANT)
    {

      if( t == 1)
        {
          n_yhat[1]=normal_rng(last_y[index+1]*beta+alpha , sigma);
          
        }
      else
        {
          n_yhat[t]=normal_rng(n_y[t-1]*beta+alpha , sigma);
        }
    }
  
}



