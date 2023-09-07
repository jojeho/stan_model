data{

  int T;
  int N;
  int R;
  array[T] matrix[N,R] Xs;

}

parameters{

  real mu ;
  real<lower=0.01> sigma;

 }


model{

  mu ~normal(0,1);
  sigma ~ normal(0,2);
  
  for(t in 1:T)
    {
      for(n in 1:N)
        {
          for(r in 1:R)
            {
              Xs[t][n,r] ~ normal(mu,sigma);
            }
        }
        
    }

}

