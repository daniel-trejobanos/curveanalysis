##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of coefficients
  int<lower=0> T; #number of timepoints
   real MAXT;
  matrix<lower=0,upper=1>[M,T] X;
  matrix[M,M] D;
  real timeN[T];
  real time[T];
  matrix[N,T] OD; #optical density data
  real MIN;
  real MAX;
  real SIGMA_MU;
  real SIGMA_A;
   real<lower=0> LAMBDA;
   int P;
   int prior[N,P];
   #real S;
}
transformed data{
  matrix[N,T] FPOD;
  for(i in 1:N)
    for(j in 1:T)
      FPOD[i,j] = log(OD[i,j]/OD[i,1]);
}
parameters {
 #real<lower=0, upper=S> s;
  matrix<lower=0,upper=1>[N,3] rate;
 real<lower=0> sigma_mu[3];
 real<lower=0> lambda[3];
 real<lower=0> MU[3];
 real <lower=0> sigma_o;
 matrix<lower=0,upper=3.1>[M,N] A_coef;

}
transformed parameters{
    matrix[P,P] lp[N]; 
    matrix[N,T] F;
    F=A_coef'*X;
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
            
            lpS_12[n,t+1]=lpS_12[n,t]+normal_lpdf(F[n,t]|rate[n,1]*MAXT*timeN[t],lambda[1]);
            lpS_21[n,t+1]=lpS_21[n,t]+normal_lpdf(F[n,t]|rate[n,2]*MAXT*timeN[t],lambda[2]);
            lpS_32[n,t+1]=lpS_32[n,t]+normal_lpdf(F[n,t]|rate[n,3]*MAXT*timeN[t],lambda[3]);
        }
        
        for(s1 in 1:P){
          log_pROW[n,s1] =  lpS_21[n,T + 1] + lpS_12[n,prior[n,s1]] - lpS_21[n,prior[n,s1]];
        }
        for(s1 in 1:P)
         for(s2 in (s1+1):P){
          lp[n,s1,s2]=  lp[n,s1,s2] + lpS_32[n,T+1] + log_pROW[n,s1]  - lpS_32[n,prior[n,s2]];
        }
      }
    }
}

model { 
  lambda ~  normal(0,LAMBDA);
  sigma_mu ~ normal(0,SIGMA_MU);
  sigma_o  ~ normal(0,0.1);
 for (n in 1:N){
      A_coef[,n]~ normal(0,SIGMA_A);
      rate[n,1]~ normal(MU[1],sigma_mu[1]);
      rate[n,2]~ normal(MU[2],sigma_mu[2]);
      rate[n,3]~ normal(MU[3],sigma_mu[3]);
 }
 
 for(n in 1:N){
   target+=log_sum_exp(to_vector(lp[n]))+normal_lpdf(FPOD[n,]|A_coef[,n]'*X,sigma_o);
 }
 
}

generated quantities{
   
  
     matrix[N,T] OD_pred;
     matrix[N,T] DOD_pred;

    for (n in 1:N){
      
      for (t in 1:T) {
      
        OD_pred[n,t] = A_coef[,n]'*X[,t];
        DOD_pred[n,t]= A_coef[,n]'*D*X[,t];
      }
  }
}
