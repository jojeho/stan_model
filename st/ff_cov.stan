
data{
  int LT;
  array[LT] real lx;
  int ST;
  array[ST] real x;

  int K;
  array[ST] vector[K] cov;
 
}

transformed data{
  int gap = LT-ST;
}

parameters{
  real l_alpha;
  real l_beta;
  real l_sigma;

  real s_alpha;
  real s_beta;
  real s_sigma;

  vector[K] cov_beta;
  
}

transformed parameters{
  array[LT-1] real ly;
  array[ST-1] real y;
  for(t in 2:LT)
    {
      ly[t-1] = l_alpha+l_beta*lx[t-1];
    }
  
  for(t in 2:ST)
    {
      y[t-1]= ly[t+gap-1] + s_alpha+s_beta*x[t-1] + dot_product(cov_beta,cov[t-1]);
    }
}

model{

  l_alpha ~ normal(0,1);
  l_beta ~ normal(0,1);
  l_sigma ~ normal(0,1);

  s_alpha ~ normal(0,1);
  s_beta ~ normal(0,1);
  s_sigma ~ normal(0,0.1);

  cov_beta ~ normal(0,2);
  
  lx[2:] ~ normal(ly,l_sigma);
  x[2:] ~ normal(y  ,s_sigma);
}


generated quantities{
  real y_pred;
  //y_pred = normal_rng(l_alpha+l_beta*lx[LT] +s_alpha+s_beta*x[ST],s_sigma);
  //y_pred = normal_rng(l_alpha+l_beta*x[LT],l_sigma);
  y_pred = s_alpha+s_beta*x[ST]+ dot_product(cov_beta,cov[ST]);
  //y_pred = normal_rng(s_alpha+s_beta*x[ST],s_sigma);
}
