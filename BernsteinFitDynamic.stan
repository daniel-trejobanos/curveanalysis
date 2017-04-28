##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of coefficients
  int<lower=0> T; #number of timepoints
   real MAXT;
  matrix<lower=0,upper=1>[M,T] X;
  matrix[M,M] D;
 # matrix[M,M] I;
  #matrix[M,M] Q;
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
   int L;
   matrix[M,L] Xs;
   matrix[M,T-L+1] Xa;
   #real S;
}
transformed data{
  real<lower=0> MI;
  matrix[N,T] FPOD;
  matrix[M,M] Identity;
   MI=M;
  for(i in 1:N)
    for(j in 1:T)
      #FPOD[i,j] = log(OD[i,j])-log(OD[i,1]);
      #FPOD[i,j] = log(OD[i,j])-MIN;
      FPOD[i,j] = OD[i,j];
      #FPOD[i,j]= OD[i,j]-OD[i,1];
  for(i in 1:M)
    for(j in 1:M)
      if(i==j)
       Identity[i,j]=1.0;
      else
       Identity[i,j]=0.0;
}
parameters {
 #real<lower=0, upper=S> s;
  matrix<lower=0, upper=MAX>[M,N] A_coef;
  real<lower=0> sigma_a[N];
  real<lower=0,upper=1> sigma_o[2];
  matrix<lower=0,upper=1>[N,2] rate;
 real<lower=0> sigma_mu[2];
 real<lower=0> lambda[2];
 real<lower=0,upper=0.5> MUG;
 real<lower=0,upper=0.1> MUE;
# cov_matrix[M] SIGMA[N];

# matrix<lower=0>[M,N] beta;
 #real<lower=0> regularize;
 
}
transformed parameters{
    matrix[P,P] lp[N]; 
    matrix[M,N] DA_coef;
    matrix[N,T] loss;
    matrix[N,T] curve;
    
     #matrix[M,N] alpha;
    DA_coef=((1/MAXT)*A_coef'*D)';
    
     curve= A_coef'*X;
     loss = (DA_coef'*X)./curve;
     #alpha[1,]=beta[1,];
     #alpha[2:M,]=beta[2:M,]-beta[1:(M-1),];
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
            
            lpS_12[n,t+1]=lpS_12[n,t]+normal_lpdf(curve[n,t]|A_coef[1,n],lambda[1]);
            lpS_21[n,t+1]=lpS_21[n,t]+normal_lpdf(loss[n,t]|rate[n,1],lambda[1]);
            lpS_32[n,t+1]=lpS_32[n,t]+normal_lpdf(loss[n,t]|rate[n,2],lambda[2]);
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
   #MU~ normal(0.25,0.25);
  lambda ~  normal(0,1);
  #sigma_mu ~ cauchy(0,SIGMA_MU);
  sigma_o~ normal(0,0.01);
  sigma_a ~ normal(0,0.1);
  sigma_mu ~ normal(0,0.1);
  #regularize ~ normal(0,1);
   
 for (n in 1:N){
      #loss[n,1:3] ~ normal(0,0.1);
      #SIGMA[n] ~ inv_wishart(M+2,0.1*Identity);
     
      rate[n,1]~ normal(MUG,sigma_mu[1]);
      rate[n,2]~ normal(MUE,sigma_mu[2]);
    #  alpha[,n] ~ normal(0.1,0.05);
      #A_coef[,n] ~ multi_normal(beta[,n],SIGMA[n]);
      A_coef[1,n] ~normal(0,0.1);
      A_coef[2:M,n]  ~ normal(0,sigma_a[n]);
 }
 
 for(n in 1:N){
   target+=log_sum_exp(to_vector(lp[n]))+normal_lpdf(FPOD[n,1:L]|A_coef[,n]'*X[,1:L],sigma_o[1])+normal_lpdf(FPOD[n,(L+1):T]|A_coef[,n]'*X[,(L+1):T],sigma_o[2]);#+normal_lpdf(FPOD[n,1:L]|A_coef[,n]'*Q*Xs,sigma_o);#+normal_lpdf(0|DA_coef[,N]'*Q,regularize);
    #target+=normal_lpdf(FPOD[n,]|A_coef[,n]'*X,sigma_o[n]);
 }
 
}

generated quantities{
   
  
     matrix[N,T] OD_pred;
     matrix[N,T] DOD_pred;

    for (n in 1:N){
      
      for (t in 1:T) {
      
        OD_pred[n,t] = A_coef[,n]'*X[,t];
        DOD_pred[n,t] = DA_coef[,n]'*X[,t];
      }
  }
}
