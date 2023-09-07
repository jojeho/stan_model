data{
  int N;
  int R;
  matrix[N,R] X;

}

parameters{

  real mu ;
  real<lower=0.01> sigma;

 }


model{

  mu ~normal(0,1);
  sigma ~ normal(0,2);
  
  for(n in 1:N)
    {
      X[n,] ~ normal(mu,sigma);
    }
}

