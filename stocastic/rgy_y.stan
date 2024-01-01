data {
  int<lower=0> T;   // # time points (equally spaced)
  vector[T] y;      // mean corrected return at time t
  int N ;
  array[N] int GT;
}
parameters {
  real mu;                     // mean log volatility
  real<lower=-1,upper=1> phi;  // persistence of volatility
  real<lower=0> sigma;         // white noise shock scale
  vector[T] h_std;
}

transformed parameters{

}
  
  
model {
  phi ~ uniform(-1, 1);
  sigma ~ cauchy(0, 5);
  h_std ~ std_normal();
  mu ~ cauchy(0, 11);
  int pos=1;
  for(n in 1:N)
    {
      int PT=GT[n];
      vector[PT] yy=segment(y,pos,PT);
      vector[PT] hh_std=segment(h_std,pos,PT);
      vector[PT]  hh = hh_std * sigma;  
      hh[1] /= sqrt(1 - phi * phi);  // rescale h[1]
      hh += mu;
      for (t in 2:PT)
        hh[t] += phi * (hh[t-1] - mu);
      for (t in 1:PT)
        yy[t] ~ normal(0, exp(hh[t] / 2));
      pos = pos + PT;
    }
}
