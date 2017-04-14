##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of coefficients
  int<lower=0> T; #number of timepoints
  matrix<lower=0,upper=1>[M,T] X;
  matrix[M,M] D;
  real timeN[T];
  matrix[N,T] OD; #optical density data
   real MU[3];
  real MIN;
  real MAX;
  real SIGMA_MU;
   
  real V;
}
transformed data{
  real<lower=0> MI;
  matrix[N,T] FPOD;
   MI=M;
  for(i in 1:N)
    for(j in 1:T)
      #FPOD[i,j] = log(OD[i,j])-log(OD[i,1]);
      FPOD[i,j] = log(OD[i,j])-MIN;
}
parameters {
 
  matrix<lower=0,upper=MAX-MIN>[M,N] A_coef;
  real<lower=0> sigma_a[N];
  real<lower=0,upper=1> sigma_o[N];
#  real intercept[N];
  matrix<lower=0>[N,3] rate;
 real<lower=0> sigma_mu;
 real<lower=0> LAMBDA;
}
transformed parameters{
    matrix[M-2,M-2] lp[N]; 
    matrix[M-1,N] DA_coef;
    
    for(n in 1:N){
      for(i in 2:M)
        DA_coef[i-1,n]= A_coef[i,n]-A_coef[i-1,n];
    lp[n] = rep_matrix(-log(M-2),M-2,M-2); 
    for (s1 in 2:M-1) 
      for (s2 in 2:M-1) 
      for (t in 2:M-1){
        if(t<s1){
          lp[n,s1-1,s2-1] = lp[n,s1-1,s2-1] + normal_lpdf(DA_coef[t-1,n]|rate[n,1],LAMBDA);
          #exponential_lpdf((DA_coef[t-1,n]-rate[n,1])^2|LAMBDA);
        }else{
          if(t<s2){
            lp[n,s1-1,s2-1] = lp[n,s1-1,s2-1] + normal_lpdf(DA_coef[t-1,n]|rate[n,2],LAMBDA);
            #exponential_lpdf((DA_coef[t-1,n]-rate[n,2])^2|LAMBDA);
          }else{
            lp[n,s1-1,s2-1] = lp[n,s1-1,s2-1] + normal_lpdf(DA_coef[t-1,n]|rate[n,3],LAMBDA);
            #exponential_lpdf((DA_coef[t-1,n]-rate[n,3])^2|LAMBDA);
          }
        }
      } 
    }
}

model { 
  LAMBDA ~  cauchy(0,10);
  sigma_mu ~ cauchy(0,SIGMA_MU);
  sigma_o~ cauchy(0,0.1);
   sigma_a ~ cauchy(0,0.1);
  rate[,1]~ normal(MU[1],sigma_mu);
   rate[,2]~ normal(MU[2],sigma_mu);
   rate[,3]~ normal(MU[3],sigma_mu);
   #rate[,1]~ cauchy(0,SIGMA_MU);
   #rate[,2]~ cauchy(0,SIGMA_MU);
   #rate[,3]~ cauchy(0,SIGMA_MU);
 for (n in 1:N){
   A_coef[,n]  ~ normal(0,sigma_a[n]);
 }
 
 for(n in 1:N){
    target+=log_sum_exp(to_vector(lp[n]))+normal_lpdf(FPOD[n,]|A_coef[,n]'*X,sigma_o[n]);
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
        OD_pred[n,t] = normal_rng(A_coef[,n]'*X[,t],sigma_o[n]);
        DOD_pred[n,t] = normal_rng(A_coef[,n]'*D*X[,t],sigma_o[n]);
      }
  }
}
