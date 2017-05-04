##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of coefficients
  int<lower=0> T; #number of timepoints
   real MAXT;
  matrix<lower=0,upper=1>[M,T] X;
  matrix[M,M] D;
  #matrix[M,M] Id;
  #matrix[M,M] Q;
  vector[T] timeN;
  real time[T];
  matrix[N,T] OD; #optical density data
  real MIN;
  real MAX;
  real SIGMA_MU;
  real SIGMA_A;
   real<lower=0> LAMBDA;
    int P;
   int prior[P];
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
  matrix<lower=0,upper=MAX>[M,N] A_coef;
  real<lower=0> sigma_a[N];
  real<lower=0,upper=1> sigma_o[2];
  vector<lower=0,upper=1>[3] rate;
 real<lower=0> sigma_mu[3];
 real<lower=0> lambda[3];
 real<lower=0> MUG;
 real<lower=0> MUE;
 real<lower=0> MUL;
# cov_matrix[M] SIGMA[N];

# matrix<lower=0>[M,N] beta;
 #real<lower=0> regularize;
 
}
transformed parameters{
    matrix[P,P] lp[N]; 
    matrix[N,T] loss;
    matrix[N,T] curve;
    matrix[N,T] derivative;
        
    matrix[M,N] DA_coef;
   # real<lower=0> beta;
     #matrix[M,N] alpha;
   
    DA_coef=((1/MAXT)*A_coef'*D)';
     curve=A_coef'*X;
     derivative=DA_coef'*X;
     loss = derivative ./curve;
    # beta= MUG-MUE;
     #alpha[1,]=beta[1,];
     #alpha[2:M,]=beta[2:M,]-beta[1:(M-1),];
    {
      row_vector[T] logCurve;
      row_vector[T] g1;
      real b;
      real b0;
      for(n in 1:N){
        lp[n] = rep_matrix(-log(P),P,P);
        logCurve=log(curve[n,]/curve[n,1]);
        for(s1 in 1:4){
         for(s2 in (s1+1):P){
           g1=(timeN-rep_vector(timeN[prior[s2]],T))';
           b0=logCurve[1];
           b=logCurve[prior[s2]];
           for(i in 1:T)
              if(i<prior[s1])
                lp[n,s1,s2]=  lp[n,s1,s2] + normal_lpdf(logCurve[i]|MAXT*rate[1]*timeN[i]+b0,lambda[1]) ;
            else{
               
              if(i<prior[s2]){
                lp[n,s1,s2]=  lp[n,s1,s2] +
          normal_lpdf(logCurve[i] |MAXT*rate[2]*g1[i]+b,lambda[2]);
              }else{
                lp[n,s1,s2]=  lp[n,s1,s2] + normal_lpdf(logCurve[i]|MAXT*rate[3]*g1[i]+b,lambda[3]);
              }
                
              }
               
            }
          
        }
      }
    }
    
}

model { 
  lambda ~  normal(0,LAMBDA);
  sigma_o~ normal(0,0.01);
  sigma_a ~ normal(0,0.1);
  sigma_mu ~ normal(0,SIGMA_MU);
   rate[1]~ normal(MUL,sigma_mu[1]);
    rate[2]~ normal(MUG,sigma_mu[2]);
    rate[3]~ normal(MUE,sigma_mu[3]);
 for (n in 1:N){
      
    
      A_coef[1,n] ~normal(0.07,0.01);
      A_coef[2:M,n]  ~ normal(0,sigma_a[n]);
 }
 
 for(n in 1:N){
   target+=log_sum_exp(to_vector(lp[n]))+normal_lpdf(FPOD[n,1:L]|A_coef[,n]'*X[,1:L],sigma_o[1])+normal_lpdf(FPOD[n,(L+1):T]|A_coef[,n]'*X[,(L+1):T],sigma_o[2]);
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
