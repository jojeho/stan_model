
  int pos=1;
  for(i in 1:GI)
    {
      int T=GT[i];
      matrix[K,T] ob=rep_matrix(0,K,T);
      for(t in 1:T)
        {
          for(k in 1:K)
            {
              for(n in 1:N)
                {
                  ob[k,t] +=normal_lpdf(X[n,t+pos-1]|mu[k] ,sigma);
                }
            }
        }
      
      target += hmm_marginal(ob,Gamma,rho);
      pos += T;
    }