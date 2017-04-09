##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of basis functions 
  int<lower=0> T; #number of timepoints
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
   matrix<lower=0,upper=1>[N,T] rate;
 real<lower=0> lambda;
}
transformed parameters{
   real<lower=0> tau[N];
   for(i in 1:N)
    tau[i]=trace((A_coef[,i]'*(D*X)-rate[i,] .* (A_coef[,i]'*X))'*(A_coef[,i]'*(D*X)-rate[i,] .* (A_coef[,i]'*X)));
   #tau= trace((A_coef'*(D*X)-rate .* (A_coef'*X))'*(A_coef'*(D*X)-rate .* (A_coef'*X)));
    
}

model {
  
 lambda~cauchy(0,5);
 for (n in 1:N){
   sigma_a[n]~cauchy(0,0.1);
  sigma_o[n]~ cauchy(0,0.1);
   tau[n]~exponential(lambda);
   A_coef[,n]~normal(0,sigma_a[n]);
    for (t in 1:T) {
      
       OD[n,t] ~ normal(A_coef[,n]'*to_vector(X[,t]),sigma_o[n]);
       
    }
     
 }
}

generated quantities{
   matrix[N,T] OD_pred;
    for (n in 1:N)
   for (t in 1:T) {
      
       OD_pred[n,t] = normal_rng(A_coef[,n]'*to_vector(X[,t]),sigma_o[n]);
       
    }
}