##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of basis functions 
  int<lower=0> T; #number of timepoints
  #int<lower=1> K; 
  #real lambda;  
  matrix[N,T] OD; #optical density data
   matrix[M,T] X; #basis functions
  matrix[M,M] D;
  real<lower=0,upper=1> MU[2];
  real LAMBDA;
  real SIGMA_MU;
  real V;
}
transformed data{
  matrix[N,T] logOD;
  for(i in 1:N)
    for(j in 1:T)
      logOD[i,j] = log(OD[i,j]);
}
parameters {
  matrix<lower=0,upper=1>[M,N] A_coef;
  #row_vector<lower=0,upper=1>[T] rate;
  real<lower=0> sigma_a[N];
  real<lower=0,upper=1> sigma_o[N];
   matrix<lower=0,upper=1>[N,3] rate;
 real<lower=0> lambda;
 #simplex[K] pi;
 #real<lower=0,upper=1> mu[2];
 real<lower=0> sigma_mu[2];
 
}
transformed parameters{
    matrix[M,M] lp[N]; 
  
    
    for(n in 1:N){
    lp[n] = rep_matrix(-log(M),M,M); 
    for (s1 in 1:M) 
      for (s2 in 1:M) 
      for (t in 1:M){
        if(t<s1){
          lp[n,s1,s2] <-lp[n,s1,s2] + exponential_log(((A_coef[,n]'*D[,t])- rate[n,1] * A_coef[t,n])^2, LAMBDA)+normal_lpdf(0|(A_coef[,n]'*D[,t]),V);
        }else{
          if(t<s2){
            lp[n,s1,s2] <-lp[n,s1,s2] + exponential_log(((A_coef[,n]'*D[,t])- rate[n,2] * A_coef[t,n])^2, LAMBDA);
          }else{
            lp[n,s1,s2] <-lp[n,s1,s2] +exponential_log(((A_coef[,n]'*D[,t])- rate[n,3] * A_coef[t,n])^2, LAMBDA);
          }
        }
      } 
    }
}

model {

 #lambda~cauchy(0,LAMBDA);
 for (n in 1:N){
   sigma_a[n]~normal(0,0.1);
  sigma_o[n]~ normal(0,0.1);
   A_coef[,n]~normal(0,sigma_a[n]);
   
   for(k in 1:2){
      sigma_mu[k] ~ normal(0,SIGMA_MU);
      rate[n,k]~ normal(MU[k],sigma_mu[k]);
      #rate[n,k]~ normal(mu[k],0.01);
      
      
   }
    target+=log_sum_exp(to_vector(lp[n])) + normal_lpdf(OD[n,]|A_coef[,n]'*X,sigma_o[n]);
 }
}

generated quantities{
   matrix[N,T] OD_pred;
   matrix[N,T] rate_pred;
   #int<lower=1,upper=M> tau[N];
    #simplex[M] sp;
    for (n in 1:N){
   for (t in 1:T) {
      
       OD_pred[n,t] = normal_rng(A_coef[,n]'*to_vector(X[,t]),sigma_o[n]);
       rate_pred[n,t]= (A_coef[,n]'*(D*to_vector(X[,t])))/(A_coef[,n]'*to_vector(X[,t]));
    }
    
   # sp = softmax(lp[,n]);
    #tau[n] = categorical_rng(sp);
  }
}
