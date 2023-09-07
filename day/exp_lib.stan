  vector gp_pred_rng(real[] x2,
                     vector y1,
                     real[] x1,
                     real alpha,
                     real rho,
                     real sigma,
                     real delta) {
    int N1 = rows(y1);
    int N2 = size(x2);
    vector[N2] f2;
    {
      matrix[N1, N1] L_K;
      vector[N1] K_div_y1;
      matrix[N1, N2] k_x1_x2;
      matrix[N1, N2] v_pred;
      vector[N2] f2_mu;
      matrix[N2, N2] cov_f2;
      matrix[N2, N2] diag_delta;
      matrix[N1, N1] K;
      K = cov_exp_quad(x1, alpha, rho);
      for (n in 1:N1)
        K[n, n] = K[n,n] + square(sigma);
      L_K = cholesky_decompose(K);
      K_div_y1 = mdivide_left_tri_low(L_K, y1);
      K_div_y1 = mdivide_right_tri_low(K_div_y1', L_K)';
      k_x1_x2 = cov_exp_quad(x1, x2, alpha, rho);
      f2_mu = (k_x1_x2' * K_div_y1);
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);
      cov_f2 = cov_exp_quad(x2, alpha, rho) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(delta, N2));

      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }

matrix gp(array[] real x , real tau , real smooth,real sigma)
  {
    int T=size(x);
    matrix[T, T] L_K;
    matrix[T, T] K = cov_exp_quad(x, tau, smooth);
    real sq_sigma = square(sigma);

    // diagonal elements
    for (n1 in 1:T)
      K[n1, n1] = K[n1, n1] + sq_sigma;

    L_K = cholesky_decompose(K);
    return L_K;
  }


matrix gp_sin(array[] real x , real tau , real smooth,real sigma)
  {
    int T=size(x);
    matrix[T, T] L_K;
    matrix[T, T] K = cov_exp_quad(x, tau, smooth);
    real sq_sigma = square(sigma);

    // diagonal elements
    for (n1 in 1:T)
      K[n1, n1] = K[n1, n1] + sq_sigma;

    L_K = cholesky_decompose(K);
    return L_K;
  }


  // void gp_lp(real tau_prior,real smooth,real sigma)
  // {
  //   tau ~ normal(0,tau_prior);
  //   smooth~normal(0,smooth_prior);
  //   sigma ~std_normal();
  // }
