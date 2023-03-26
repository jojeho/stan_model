functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]); 
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) / 
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) / 
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines 
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) + 
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }

  vector spline(real[] X ,row_vector a_raw,real a0 , real tau,matrix B,int num_basis)
  {
    int num_data=size(X);
    row_vector[num_basis] a; // spline coefficients
    vector[num_data] Y_hat;
    a = a_raw*tau;
    Y_hat = a0*to_vector(X) + to_vector(a*B);
    return Y_hat;
  }

  matrix basis(real [] X ,int num_data,int num_knots,
               vector knots,int spline_degree)
  {
    int num_basis = num_knots + spline_degree - 1; // total number of B-splines

    matrix[num_basis, num_data] B;  // matrix of B-splines
    vector[spline_degree + num_knots] ext_knots_temp;
    vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
    ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
    ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
    for (ind in 1:num_basis)
      B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
    B[num_knots + spline_degree - 1, num_data] = 1;

    return B;
  }

  vector[] spline_oblik_k(vector[] Y,real[] x ,int K,vector lambda,row_vector[] a_raw , vector a0 , vector tau
                          ,matrix B,int num_basis,int num_data,
                          real sigma
                          )
  {
    int T = size(Y);
    vector[K] alpha[T];

    for(t in 1:T)
    {
      alpha[t,]=log(lambda);
      for(k in 1:K)
        {
          vector[num_data] y_hat;
          y_hat = spline(x,a_raw[k],a0[k],tau[k],B,num_basis);
          alpha[t,k] += normal_lpdf(Y[t]|y_hat, sigma);
        }
    }

    return alpha;
  }

  void spline_lp(row_vector a_raw,real a0,real tau,real sigma)

  {
      a_raw ~ normal(0, 5);
      a0 ~ cauchy(0, 5);
      tau ~ cauchy(0, 5);
      sigma ~ normal(0,5);
  }
}


data {
  int num_data;             // number of data points
  int num_knots;            // num of knots
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)
  real Y[num_data];
  real X[num_data];

  
}

transformed data {
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, num_data] B;  // matrix of B-splines
  B = basis(X,num_data,num_knots,knots,spline_degree);
}

parameters {
  row_vector[num_basis] a_raw; 
  real a0;  // intercept
  real<lower=0> sigma; 
  real<lower=0> tau;   
}

transformed parameters {

  vector[num_data] Y_hat;
  Y_hat=spline(X,a_raw,a0,tau,B,num_basis);
}

model {
  // Priors
  spline_lp(a_raw ,a0,sigma,tau);
  //Likelihood
  Y ~ normal(Y_hat, sigma);
}


