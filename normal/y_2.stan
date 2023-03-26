functions{
  void normal_lp(real mu)
  {
    mu ~ normal(0,10);
  }

  real v_lpdf(vector x , real mu,real sigma)
  {
    return normal_lpdf(x|mu,sigma);
  }
    
}

data{
  int N;
  vector[N] x;
  int v;
}

parameters{
  real mu;
  real<lower=0> sigma;
}

model{

  mu ~ normal(0,1);
  sigma ~normal(0,1);
  x ~ student_t(v,mu,sigma);
  //x ~ v(mu,sigma);
}

generated quantities{
  array[N] real oblik;
  array[N] real y_hat;
  for(t in 1:N)
    {
      oblik[t]=student_t_lpdf(x|v,mu,sigma);
      y_hat[t]=student_t_rng(v,mu,sigma);
    }
}
