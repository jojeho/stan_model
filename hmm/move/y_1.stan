data {
  int<lower=0> T; // length of the time series
  int ID[T]; // track identifier
  vector[T] steps; // step lengths
  vector[T] angles; // turning angles
  int<lower=1> N; // number of states
  int nCovs; // number of covariates
  matrix[T,nCovs+1] covs; // covariates
}


parameters {
  positive_ordered[N] mu; // mean of gamma - ordered
  vector<lower=0>[N] sigma; // SD of gamma
  // unconstrained angle parameters
  vector[N] xangle;
  vector[N] yangle;
  // regression coefficients for transition probabilities
  matrix[N*(N-1),nCovs+1] beta;
}


transformed parameters {
  vector<lower=0>[N] shape;
  vector<lower=0>[N] rate;
  vector<lower=-pi(),upper=pi()>[N] loc;
  vector<lower=0>[N] kappa;
  // derive turning angle mean and concentration
  for(n in 1:N) {
    loc[n] = atan2(yangle[n], xangle[n]);
    kappa[n] = sqrt(xangle[n]*xangle[n] + yangle[n]*yangle[n]);
  }
  // transform mean and SD to shape and rate
  for(n in 1:N)
    shape[n] = mu[n]*mu[n]/(sigma[n]*sigma[n]);
  for(n in 1:N)
    rate[n] = mu[n]/(sigma[n]*sigma[n]);
}
model {
  vector[N] logp;
  vector[N] logptemp;
  matrix[N,N] gamma[T];
  matrix[N,N] log_gamma[T];
  matrix[N,N] log_gamma_tr[T];
  // priors
  mu ~ normal(0, 5);
  sigma ~ student_t(3, 0, 1);
  xangle[1] ~ normal(-0.5, 1); // equiv to concentration when yangle = 0
  xangle[2] ~ normal(2, 2);
  yangle ~ normal(0, 0.5); // zero if mean angle is 0 or pi
  // derive array of (log-)transition probabilities
  for(t in 1:T) {
    int betarow = 1;
    for(i in 1:N) {
      for(j in 1:N) {
        if(i==j) {
          gamma[t,i,j] = 1;
        } else {
          gamma[t,i,j] = exp(beta[betarow] * to_vector(covs[t]));
          betarow = betarow + 1;
        }

      }
    }
    // each row must sum to 1
    for(i in 1:N)
      log_gamma[t][i] = log(gamma[t][i]/sum(gamma[t][i]));
  }
  // transpose
  for(t in 1:T)
    for(i in 1:N)
      for(j in 1:N)
        log_gamma_tr[t,j,i] = log_gamma[t,i,j];
  // likelihood computation
  for (t in 1:T) {
    // initialise forward variable if first obs of track
    if(t==1 || ID[t]!=ID[t-1])
      logp = rep_vector(-log(N), N);
    for (n in 1:N) {
      logptemp[n] = log_sum_exp(to_vector(log_gamma_tr[t,n]) + logp);
      if(steps[t]>=0)
        logptemp[n] = logptemp[n] + gamma_lpdf(steps[t] | shape[n], rate[n]);
      if(angles[t]>=(-pi()))
        logptemp[n] = logptemp[n] + von_mises_lpdf(angles[t] | loc[n], kappa[n]);
    }
    logp = logptemp;
    // add log forward variable to target at the end of each track
    if(t==T || ID[t+1]!=ID[t])
      target += log_sum_exp(logp);
  }
}


generated quantities {
  int<lower=1,upper=N> viterbi[T];
  real stateProbs[T,N];
  vector[N] lp;
  vector[N] lp_p1;
  // Viterbi algorithm (most likely state sequence)
  {
    real max_logp;
    int back_ptr[T, N];
    real best_logp[T, N];
    for (t in 1:T) {
      if(t==1 || ID[t]!=ID[t-1]) {
        for(n in 1:N)
          best_logp[t, n] = gamma_lpdf(steps[t] | shape[n], rate[n]);
      } else {
        for (n in 1:N) {
          best_logp[t, n] = negative_infinity();
          for (j in 1:N) {
            real logp;
            logp = best_logp[t-1, j] + log_theta[t,j,n];
            if(steps[t]>0)
              logp = logp + gamma_lpdf(steps[t] | shape[n], rate[n]);
            if(angles[t]>(-pi()))
              logp = logp + von_mises_lpdf(angles[t] | loc[n], kappa[n]);
            if (logp > best_logp[t, n]) {
              back_ptr[t, n] = j;
              best_logp[t, n] = logp;
            }
          }
        }
      }
    }
    for(t0 in 1:T) {
      int t = T - t0 + 1;
      if(t==T || ID[t+1]!=ID[t]) {
        max_logp = max(best_logp[t]);

 
        for (n in 1:N)
          if (best_logp[t, n] == max_logp)
            viterbi[t] = n;
      } else {
        viterbi[t] = back_ptr[t+1, viterbi[t+1]];
      }
    }
  }
  // forward-backward algorithm (state probabilities)
  {
    real logalpha[T,N];
    real logbeta[T,N];
    real llk;
    // log alpha probabilities
    for(t in 1:T) {
      if(t==1 || ID[t]!=ID[t-1]) {
        for(n in 1:N)
          lp[n] = -log(N);
      }
      for (n in 1:N) {
        lp_p1[n] = log_sum_exp(to_vector(log_theta_tr[t,n]) + lp);
        if(steps[t]>=0)
          lp_p1[n] = lp_p1[n] + gamma_lpdf(steps[t] | shape[n], rate[n]);
        if(angles[t]>=(-pi())) {
          lp_p1[n] = lp_p1[n] + von_mises_lpdf(angles[t] | loc[n], kappa[n]);
        }
        logalpha[t,n] = lp_p1[n];
      }
      lp = lp_p1;
    }
    // log beta probabilities
    for(t0 in 1:T) {
      int t = T - t0 + 1;
      if(t==T || ID[t+1]!=ID[t]) {
        for(n in 1:N)
          lp_p1[n] = 0;
      } else {
        for(n in 1:N) {
          lp_p1[n] = log_sum_exp(to_vector(log_theta_tr[t+1,n]) + lp);
          if(steps[t+1]>=0)
            lp_p1[n] = lp_p1[n] + gamma_lpdf(steps[t+1] | shape[n], rate[n]);
          if(angles[t+1]>=(-pi()))
            lp_p1[n] = lp_p1[n] + von_mises_lpdf(angles[t+1] | loc[n], kappa[n]);
        }
      }
      lp = lp_p1;
      for(n in 1:N)
        logbeta[t,n] = lp[n];
    }for (n in 1:N)
       if (best_logp[t, n] == max_logp)
         viterbi[t] = n;
  } else {
    viterbi[t] = back_ptr[t+1, viterbi[t+1]];
  }
}
}
// forward-backward algorithm (state probabilities)
{
  real logalpha[T,N];
  real logbeta[T,N];
  real llk;
  // log alpha probabilities
  for(t in 1:T) {
    if(t==1 || ID[t]!=ID[t-1]) {
      for(n in 1:N)
        lp[n] = -log(N);
    }
    for (n in 1:N) {
      lp_p1[n] = log_sum_exp(to_vector(log_theta_tr[t,n]) + lp);
      if(steps[t]>=0)
        lp_p1[n] = lp_p1[n] + gamma_lpdf(steps[t] | shape[n], rate[n]);
      if(angles[t]>=(-pi())) {
        lp_p1[n] = lp_p1[n] + von_mises_lpdf(angles[t] | loc[n], kappa[n]);
      }
      logalpha[t,n] = lp_p1[n];
    }
    lp = lp_p1;
  }
  // log beta probabilities
  for(t0 in 1:T) {
    int t = T - t0 + 1;
    if(t==T || ID[t+1]!=ID[t]) {
      for(n in 1:N)
        lp_p1[n] = 0;
    } else {
      for(n in 1:N) {
        lp_p1[n] = log_sum_exp(to_vector(log_theta_tr[t+1,n]) + lp);
        if(steps[t+1]>=0)
          lp_p1[n] = lp_p1[n] + gamma_lpdf(steps[t+1] | shape[n], rate[n]);
        if(angles[t+1]>=(-pi()))
          lp_p1[n] = lp_p1[n] + von_mises_lpdf(angles[t+1] | loc[n], kappa[n]);
      }
    }
    lp = lp_p1;
    for(n in 1:N)
      logbeta[t,n] = lp[n];
  }

  for(t0 in 1:T) {
    int t = T - t0 + 1;
    if(t==T || ID[t+1]!=ID[t])
      llk = log_sum_exp(logalpha[t]);
    for(n in 1:N)
      stateProbs[t,n] = exp(logalpha[t,n] + logbeta[t,n] - llk);
  }
}
}
