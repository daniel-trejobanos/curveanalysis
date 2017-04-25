##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of coefficients
  int<lower=0> T; #number of timepoints
   real MAXT;
  matrix<lower=0,upper=1>[M,T] X;
  
  real timeN[T];
  real time[T];
  matrix[N,T] OD; #optical density data
  real MIN;
  real MAX;
  real SIGMA_MU;
   real<lower=0> LAMBDA;
   int P;
   int prior[N,P];
   #real S;
}
transformed data{
  real<lower=0> MI;
  matrix[N,T] FPOD;
  for(i in 1:N)
    for(j in 1:T)
      FPOD[i,j] = log(OD[i,j]/OD[i,1]);
}
parameters {
 #real<lower=0, upper=S> s;
  matrix<lower=0, upper=MAX>[M,N] A_coef;
  matrix<lower=0,upper=1>[N,3] rate;
 real<lower=0> sigma_mu[3];
 real<lower=0> lambda[3];
 real<lower=0> MU[3];

}
transformed parameters{
    matrix[P,P] lp[N]; 

    {
      matrix[N,T+1] lpS_12;
      matrix[N,T+1] lpS_21;
      matrix[N,T+1] lpS_32;
      matrix[N,P] log_pROW;
      real temp;
      for(n in 1:N){
        lp[n] = rep_matrix(-log(P),P,P);
        lpS_12[n,1]=0;
        lpS_21[n,1]=0;
        lpS_32[n,1]=0;
        for (t in 1:T){
            
            lpS_12[n,t+1]=lpS_12[n,t]+normal_lpdf(FPOD[,t]|rate[n,1]*timeN[t],lambda[1]);
            lpS_21[n,t+1]=lpS_21[n,t]+normal_lpdf(FPOD[,t]|rate[n,2]*timeN[t],lambda[2]);
            lpS_32[n,t+1]=lpS_32[n,t]+normal_lpdf(FPOD[,t]|rate[n,3]*timeN[t],lambda[3]);
        }
        
        for(s1 in 1:P){
          log_pROW[n,s1] =  lpS_21[n,T + 1] + lpS_12[n,prior[n,s1]] - lpS_21[n,prior[n,s1]];
        }
        for(s1 in 1:P)
         for(s2 in 2:P){
          lp[n,s1,s2]=  lp[n,s1,s2] + lpS_32[n,T+1] + log_pROW[n,s1]  - lpS_32[n,prior[n,s2]];
        }
      }
    }
}

model { 
  lambda ~  normal(0,LAMBDA);
  sigma_mu ~ cauchy(0,SIGMA_MU);
 for (n in 1:N){
      rate[n,1]~ normal(MU[1],sigma_mu[1]);
      rate[n,2]~ normal(MU[2],sigma_mu[2]);
      rate[n,3]~ normal(MU[3],sigma_mu[3]);
 }
 
 for(n in 1:N){
   target+=log_sum_exp(to_vector(lp[n]));
 }
 
}

generated quantities{
   
  
     matrix[N,T] OD_pred;
     matrix[N,T] DOD_pred;

    for (n in 1:N){
      
      for (t in 1:T) {
      
        OD_pred[n,t] = rate[n,1]'*timeN[t];
        DOD_pred[n,t]= rate[n,2]'*timeN[t];
      }
  }
}
