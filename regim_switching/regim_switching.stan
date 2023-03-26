functions{
  matrix make_oblik( vector y
                   ,vector alpha ,real rho
                   ,real y_tm1_init,vector sigma)
  {
    int T = size(y);
    matrix[T, 2] eta;
    for(t in 1:T) {
      eta[t,1] = exp(normal_lpdf(y[t]| alpha[1], sigma[1]));
      if(t==1) {
        eta[t,2] = exp(normal_lpdf(y[t]| alpha[2] + rho * y_tm1_init, sigma[2]));
      } else {
        eta[t,2] = exp(normal_lpdf(y[t]| alpha[2] + rho * y[t-1], sigma[2]));
      }
    }

    return eta;

  }

  real fill_prob(matrix eta
                 ,real xi1_init,vector p
                 )
  {
    int T = dims(eta)[1];
    vector[T] f;
    matrix[T, 2] xi;
    for(t in 1:T) {
      // for the first observation
      if(t==1) {
        f[t] = p[1]*xi1_init*eta[t,1] + // stay in state 1
          (1 - p[1])*xi1_init*eta[t,2] + // transition from 1 to 2
          p[2]*(1 - xi1_init)*eta[t,2] + // stay in state 2 
          (1 - p[2])*(1 - xi1_init)*eta[t,1]; // transition from 2 to 1
      
        xi[t,1] = (p[1]*xi1_init*eta[t,1] +(1 - p[2])*(1 - xi1_init)*eta[t,1])/f[t];
        xi[t,2] = 1.0 - xi[t,1];
    
      } else {
        //  and for the rest
      
        f[t] = p[1]*xi[t-1,1]*eta[t,1] + // stay in state 1
          (1 - p[1])*xi[t-1,1]*eta[t,2] + // transition from 1 to 2
          p[2]*xi[t-1,2]*eta[t,2] + // stay in state 2 
          (1 - p[2])*xi[t-1,2]*eta[t,1]; // transition from 2 to 1
      
        // work out xi
      
        xi[t,1] = (p[1]*xi[t-1,1]*eta[t,1] +(1 - p[2])*xi[t-1,2]*eta[t,1])/f[t];
      
        // there are only two states so the probability of the other state is 1 - prob of the first
        xi[t,2] = 1.0 - xi[t,1];
      }
    }

    return sum(log(f));
  }

  matrix fill_xi(matrix eta
                 ,real xi1_init,vector p
                 )
  {
    int T = dims(eta)[1];
    vector[T] f;
    matrix[T, 2] xi;
    for(t in 1:T) {
      // for the first observation
      if(t==1) {
        f[t] = p[1]*xi1_init*eta[t,1] + // stay in state 1
          (1 - p[1])*xi1_init*eta[t,2] + // transition from 1 to 2
          p[2]*(1 - xi1_init)*eta[t,2] + // stay in state 2 
          (1 - p[2])*(1 - xi1_init)*eta[t,1]; // transition from 2 to 1
      
        xi[t,1] = (p[1]*xi1_init*eta[t,1] +(1 - p[2])*(1 - xi1_init)*eta[t,1])/f[t];
        xi[t,2] = 1.0 - xi[t,1];
    
      } else {
        //  and for the rest
      
        f[t] = p[1]*xi[t-1,1]*eta[t,1] + // stay in state 1
          (1 - p[1])*xi[t-1,1]*eta[t,2] + // transition from 1 to 2
          p[2]*xi[t-1,2]*eta[t,2] + // stay in state 2 
          (1 - p[2])*xi[t-1,2]*eta[t,1]; // transition from 2 to 1
      
        // work out xi
      
        xi[t,1] = (p[1]*xi[t-1,1]*eta[t,1] +(1 - p[2])*xi[t-1,2]*eta[t,1])/f[t];
      
        // there are only two states so the probability of the other state is 1 - prob of the first
        xi[t,2] = 1.0 - xi[t,1];
      }
    }

    return xi;
  }

}

// saved as regime_switching_model.stan
data {
  int T;
  int N;
  array[T] vector[N] y;
}
parameters {

  array[T] vector<lower = 0, upper = 1>[2] p;
  array[T] real<lower = 0> rho;
  array[T] vector[2] alpha;
  array[T] vector<lower = 0>[2] sigma;
  array[T] real<lower = 0, upper = 1> xi1_init; 
  array[T] real y_tm1_init;
  
}
transformed parameters {

  // fill in etas
  // work out likelihood contributions
  
}
model {
  // priors
  for(t in 1:T)
    {
      p[t] ~ beta(10, 2);
      rho[t] ~ normal(0, 5);
      alpha[t] ~ normal(0, 0.05);
      sigma[t] ~ cauchy(0, 0.1);
      xi1_init[t] ~ beta(2, 2);
      y_tm1_init[t] ~ normal(0, .1);
    }
  

  for(t in 1:T)
    {
      matrix[N,2] eta;
      eta=make_oblik(y[t],alpha[t],rho[t],y_tm1_init[t],sigma[t]);
      target +=fill_prob(eta,xi1_init[t],p[t]);
    }
}


/* generated quantities{ */
/*   array[T] real y_hat; */
/*   matrix[N, 2] xi; */
/*   real f1; */
/*   real f2; */
/*   vector[2] pp; */
  
/*   for(t in 1:T) */
/*     { */
/*       matrix[N,2] eta=fill_eta(y[t],alpha[t],rho[t],y_tm1_init[t],sigma[t]); */
/*       xi = fill_xi(eta,xi1_init[t],p[t]); */
/*       pp=p[t]; */
/*       f1 = pp[1]*xi[N,1] +(1- pp[2])*xi[N,2]; */
/*       f2 = (1-pp[1])*xi[N,2] +(pp[2])*xi[N,2]; */

/*       f1 = f1/(f1+f2); */
/*       f2 = f2/(f1+f2); */

/*       y_hat[t] = (alpha[t,2] + rho[t] * y[t,N])*f2 + (alpha[t,1])*f1; */
/*     } */

/* } */

  
