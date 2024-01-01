data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
}
parameters {
  real mu;                     // mean log volatility
  real<lower=-1,upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  vector[T] h_std;
  real<lower=0> v;
}

transformed parameters{
 vector[T] h = h_std * sigma;  // now h ~ normal(0, sigma)
  h[1] /= sqrt(1 - phi * phi);  // rescale h[1]
  h += mu;
  for (t in 2:T)
    h[t] += phi * (h[t-1] - mu);
}
  
  
model {
  phi ~ uniform(-1, 1);
  sigma ~ cauchy(0, 5);
  h_std ~ std_normal();
  mu ~ cauchy(0, 11);
  for (t in 1:T)
    y[t] ~ normal(0, exp(h[t] / 2));
}
