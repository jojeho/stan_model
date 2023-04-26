generated quantities{
  matrix[K,D] prob;

  {
    int pos=1;
    for(i in 1:GI)
      {
        int T=GT[i];
        matrix[K,T] ob=rep_matrix(0,K,T);
        matrix[K,T] local_prob=rep_matrix(0,K,T);
        for(t in 1:T)
          {
            for(k in 1:K)
              {
                for(n in 1:N)
                  {
                    ob[k,t] +=normal_lpdf(X[n,pos+t-1]|mu[k] ,sigma);
                  }
              }
          }

        local_prob=hmm_forward_prob(ob,Gamma,rho);
        for(t in 1:GT[i])
          {
            prob[:,pos+t-1]= local_prob[:,t];
          }
        pos += GT[i];
      }
  }
  
}
