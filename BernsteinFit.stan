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
  row_vector<lower=0,upper=1>[T] rate;
  real<lower=0> sigma_a;
  real<lower=0> sigma_o;
  real<lower=0> sigma_tau;
}
transformed parameters{
 matrix[N,T] tau;
 for(i in 1:N)
  for(j in 1:T)
   tau[i,j]=A_coef[,i]'*D*X[,j]-rate[j];
}

model {
  
  sigma_a~cauchy(0,5);
  sigma_o~cauchy(0,5);
  sigma_tau~cauchy(0,5);
 for (n in 1:N){
   A_coef[,n]~normal(0,sigma_a);
    for (t in 1:T) {
       logOD[n,t] ~ normal(A_coef[,n]'*to_vector(X[,t]),sigma_o);
       tau[n,t]~normal(0,sigma_tau);
    }
     
 }
}
