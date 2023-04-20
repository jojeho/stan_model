data{
  #include "../../input/y.stan"
  real v;
  real sigma;
}


transformed data {
  real log_unif;
  log_unif = -log(T);
}
parameters {
  real mu1;
  //  real mu2;
}

transformed parameters {
  real mu2=mu1*-1;
  
  vector[T] lp;
  lp = rep_vector(log_unif, T);
  for (s in 1:T) {
    for (t in 1:T) {
      lp[s] = lp[s] + student_t_lpdf(y[t] |v, t < s ? mu1 : mu2 , sigma );
    }
  }
}


model {
  mu1 ~ normal(0,0.05);
  //mu2 ~ normal(0,0.05);
  target += log_sum_exp(lp);
}

generated quantities {
  int<lower=1, upper=T> s;
  s = categorical_logit_rng(lp);
  /* array[T] real log_lik; */
  /* for( t in 1:T) */
  /*   { */
  /*     if( t < s) */
  /*       { */
  /*         log_lik[t]=normal_lpdf(y[t] | mu1,sigma); */
  /*       } */
  /*     else */
  /*       { */
  /*         log_lik[t]=normal_lpdf(y[t] | 0,sigma); */
  /*       } */
  /*   } */
}
