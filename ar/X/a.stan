data {
  int N;
  int T;
  array[N] vector[T] X;

  int NT;
  vector[T] n_x;
}

parameters{
  array[N] real beta;
  array[N] real alpha;
  vector<lower=0>[N] sigma;
}

model{

  beta ~ normal(0,1);
  alpha  ~normal(0,1);
  sigma ~ normal(0,1);
  for(n in 1:N)
    {
      for(t in 2:T)
        X[n,t] ~ normal(beta[n]*X[n,t-1]+alpha[n],sigma[n]);
    }
}

generated quantities{
    
}

