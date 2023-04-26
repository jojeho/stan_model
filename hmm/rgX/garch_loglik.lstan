
  int pos=1;
  vector[N] sigma;

  for(i in 1:GI)
    {
      int T=GT[i];
      matrix[K,T] ob=rep_matrix(0,K,T);
      for(t in 1:T)
        {
          for(k in 1:K)
            {
              sigma[1] = sigma1;
              for (n in 2:N)
                {
                  sigma[n] = sqrt(alpha0
                                  + alpha1 * pow(X[n-1,t+pos-1] - mu[k], 2)
                                  + beta1 * pow(sigma[n-1], 2));
                }
              
              for(n in 1:N)
                {
                  ob[k,t] +=normal_lpdf(X[n,t+pos-1]|mu[k] ,sigma);
                }
              ob[k,t] -=log(N);
            }
        }
      
      target += hmm_marginal(ob,Gamma,rho);
      pos += T;
    }
