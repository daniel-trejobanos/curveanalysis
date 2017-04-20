##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of coefficients
  int<lower=0> T; #number of timepoints
   real MAXT;
  matrix<lower=0,upper=1>[M,T] X;
  matrix[M,M] D;
  real timeN[T];
  matrix[N,T] OD; #optical density data
   real MU[3];
  real MIN;
  real MAX;
  real SIGMA_MU;
   real<lower=0> LAMBDA;
    int P;
   int prior[P];
  
   #real S;
}
transformed data{
  real<lower=0> MI;
  matrix[N,T] FPOD;
   MI=M;
  for(i in 1:N)
    for(j in 1:T)
      #FPOD[i,j] = log(OD[i,j])-log(OD[i,1]);
      #FPOD[i,j] = log(OD[i,j])-MIN;
      FPOD[i,j] = OD[i,j];
}
parameters {
 #real<lower=0, upper=S> s;
  matrix<lower=0,upper=MAX-MIN>[M,N] A_coef;
  real<lower=0> sigma_a;
  real<lower=0,upper=1> sigma_o;
#  real intercept[N];
  matrix<lower=0>[N,3] rate;
 real<lower=0> sigma_mu;
 real<lower=0> lambda;
 # matrix<lower=0,upper=MAX>[M,N] jump;
#  real<lower=0,upper=MAX> k_0;
}
transformed parameters{
    matrix[P,P] lp[N]; 
    matrix[M,N] DA_coef;
    matrix[N,T] loss;
    
    
   # matrix<lower=0>[M,N] A_coef;
    
    DA_coef=(1/MAXT)*(A_coef'*D)';
    loss = (DA_coef'*X)./(A_coef'*X);
     
    
    {
      vector[T+1] lpS_12;
      vector[T+1] lpS_21;
      vector[T+1] lpS_32;
      vector[P] log_pROW;
      for(n in 1:N){
        lp[n] = rep_matrix(-log(P),P,P);
        lpS_12[1]=0;
        lpS_21[1]=0;
        lpS_32[1]=0;
        for (t in 1:T){
            
            lpS_12[t+1]=lpS_12[t]+normal_lpdf(loss[n,t]|rate[n,1],lambda);
            lpS_21[t+1]=lpS_21[t]+normal_lpdf(loss[n,t]|rate[n,2],lambda);
            lpS_32[t+1]=lpS_32[t]+normal_lpdf(loss[n,t]|rate[n,3],lambda);
        }
        
        for(s1 in 1:P){
          log_pROW[s1] =  lpS_21[T + 1] + lpS_12[prior[s1]] - lpS_21[prior[s1]];
        }
        for(s1 in 1:P)
         for(s2 in 2:P){
          lp[n,s1,s2]=  lp[n,s1,s2] + lpS_32[T+1] + log_pROW[s1]  - lpS_32[prior[s2]];
        }
      }
    }
}

model { 
  lambda ~  normal(0,LAMBDA);
  sigma_mu ~ normal(0,SIGMA_MU);
  sigma_o~ normal(0,0.1);
  sigma_a ~ normal(0,0.1);
  rate[,1]~ normal(MU[1],sigma_mu);
   rate[,2]~ normal(MU[2],sigma_mu);
   rate[,3]~ normal(MU[3],sigma_mu);
   #rate[,1]~ cauchy(0,SIGMA_MU);
   #rate[,2]~ cauchy(0,SIGMA_MU);
   #rate[,3]~ cauchy(0,SIGMA_MU);
 for (n in 1:N){
      
#    jump[,n] ~ chi_square(k_0/M);
  
      A_coef[,n]  ~ normal(0,sigma_a);
 }
 
 for(n in 1:N){
   target+=log_sum_exp(to_vector(lp[n]))+normal_lpdf(FPOD[n,]|A_coef[,n]'*X,sigma_o);
    #target+=normal_lpdf(FPOD[n,]|A_coef[,n]'*X,sigma_o[n]);
 }
 
}

generated quantities{
   
  # matrix[N,T] rate1_pred;
   # matrix[N,T] rate2_pred;
    # matrix[N,T] rate3_pred;
     matrix[N,T] OD_pred;
     matrix[N,T] DOD_pred;
     
     
   #int<lower=1,upper=M> tau[N];
    #simplex[M] sp;
    for (n in 1:N){
      
      for (t in 1:T) {
      
     #    rate1_pred[n,t] = normal_rng(rate[n,1]*timeN[t],sigma_o[n]);
      #  rate2_pred[n,t] = normal_rng(rate[n,2]*timeN[t],sigma_o[n]);
      #  rate3_pred[n,t] =normal_rng(rate[n,3]*timeN[t],sigma_o[n]);
        OD_pred[n,t] = A_coef[,n]'*X[,t];
        DOD_pred[n,t] = DA_coef[,n]'*X[,t];
      }
  }
}
