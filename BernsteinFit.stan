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
   matrix<lower=0,upper=1>[N,2] rate;
 real<lower=0> lambda;
 #simplex[K] pi;
 real<lower=0,upper=1> mu[2];
 #real<lower=0> sigma_mu[K];
}
transformed parameters{
   matrix[M,N] log_p;
   {
   vector[M+1] log_p_e;
   vector[M+1] log_p_l;
   log_p_e[1]=0;
   log_p_l[1]=0;
   for(n in 1:N){
    for(i in 1:M){
       log_p_e[i+1] = log_p_e[i] + exponential_log(((A_coef[,n]'*D[,i])- rate[n,1] * A_coef[i,n])^2, lambda);
      log_p_l[i+1] = log_p_l[i] + exponential_log(((A_coef[,n]'*D[,i])- rate[n,2] * A_coef[i,n])^2, lambda);
   }
    log_p[,n]=(rep_vector(-log(M)+log_p_l[M+1],M)+ head(log_p_e,M) - head(log_p_l,M));    
   }    
  }
}

model {

 lambda~cauchy(0,5);
 for (n in 1:N){
   sigma_a[n]~cauchy(0,0.1);
  sigma_o[n]~ cauchy(0,0.1);
   A_coef[,n]~normal(0,sigma_a[n]);
   
   for(k in 1:2){
      #sigma_mu[k] ~ cauchy(0,0.1);
      rate[n,k]~ normal(mu[k],0.01);
      
      
   }
    target+=log_sum_exp(log_p[,n]) + normal_lpdf(OD[n,]|A_coef[,n]'*X,sigma_o[n]);
 }
}

generated quantities{
   matrix[N,T] OD_pred;
   matrix[N,T] rate_pred;
    for (n in 1:N)
   for (t in 1:T) {
      
       OD_pred[n,t] = normal_rng(A_coef[,n]'*to_vector(X[,t]),sigma_o[n]);
       rate_pred[n,t]= (A_coef[,n]'*(D*to_vector(X[,t])))/(A_coef[,n]'*to_vector(X[,t]));
    }
}