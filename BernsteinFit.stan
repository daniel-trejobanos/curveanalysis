##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of basis functions 
  int<lower=0> T; #number of timepoints
  
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
  matrix[M,N] A_coef;
  #row_vector<lower=0,upper=1>[T] rate;
  real<lower=0> sigma_a[N];
  real<lower=0,upper=1> sigma_o[N];
}
transformed parameters{
  
}

model {
  
  
 for (n in 1:N){
   sigma_a[n]~cauchy(0,0.5);
   sigma_o[n]~cauchy(0,0.5);
   A_coef[,n]~normal(0,sigma_a[n]);
    for (t in 1:T) {
       logOD[n,t] ~ normal(A_coef[,n]'*to_vector(X[,t]),sigma_o[n]);
    }
     
 }
}

generated quantities{
  matrix[N,T] rate;
   rate = A_coef'*(D*X);
}