##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of coefficients
  int<lower=0> T; #number of timepoints
   real MAXT;
  matrix<lower=0,upper=1>[M,T] X;
  matrix[M,M] D;
  vector[T] timeN;
  real time[T];
  matrix[N,T] OD; #optical density data
  real MIN;
  real MAX;
  real SIGMA_A;
  real<lower=0> LAMBDA;
  int P;
   int L;
   int prior[P];
}
transformed data{
  matrix[N,T] FPOD;
  real log_unif_P;
  log_unif_P=-log(P);
  for(i in 1:N)
    for(j in 1:T)
      FPOD[i,j] = OD[i,j];
}
parameters {
 #real<lower=0, upper=S> s;
  vector<lower=0,upper=MAX>[M] A_coef;
  real<lower=0> sigma_a;
  real<lower=0,upper=1> sigma_o[2];
 real<lower=0> lambda[3];
 real<lower=0> MUG;
 real<lower=0> MUE;
 real<lower=0> MUL;
}
transformed parameters{
    matrix[P,P] lp; 
    row_vector[T] loss;
    
    row_vector[T] derivative;
    row_vector[T] logCurve;
     row_vector[T] curve;
    vector[M] DA_coef;
    DA_coef=((1/MAXT)*A_coef'*D)';
     curve=A_coef'*X;
     derivative=DA_coef'*X;
     loss = derivative ./curve;
     logCurve=log(curve ./ curve[1]);
        lp = rep_matrix(log_unif_P,P,P);
    {
      vector[T+1] lpS_12;
      vector[T+1] lpS_21;
      vector[T+1] lpS_32;
      vector[P] log_pROW;
      
        
        
       
        lpS_12[1]=0;
        lpS_21[1]=0;
        lpS_32[1]=0;
        for (t in 1:T){
            
            lpS_12[t+1]=lpS_12[t]+normal_lpdf(loss[t]|MUL,lambda[1]);
            lpS_21[t+1]=lpS_21[t]+normal_lpdf(loss[t]|MUG,lambda[2]);
            lpS_32[t+1]=lpS_32[t]+normal_lpdf(loss[t]|MUE,lambda[3]);
        }
        
        
        for(s1 in 1:P)
         for(s2 in (s1+1):P){
          lp[s1,s2]=  lp[s1,s2] + lpS_32[T+1] + lpS_12[prior[s1]] + (lpS_21[prior[s2]]-lpS_21[prior[s1]])  - lpS_32[prior[s2]];
        }
        
      }
    
}

model { 
  lambda ~  cauchy(0,LAMBDA);
  sigma_o~ normal(0,0.01);
  sigma_a ~ normal(0,SIGMA_A);
  MUL~ normal(0,0.01);
  MUG~ normal(0.2,0.01);
  MUE ~ normal(0.06,0.01);
 for (n in 1:N){
      
     FPOD[n,1:L] ~ normal(A_coef'*X[,1:L],sigma_o[1]);
     FPOD[n,(L+1):T]~normal(A_coef'*X[,(L+1):T],sigma_o[2]);
      A_coef[1] ~normal(0.07,0.01);
      A_coef[2:M]  ~ normal(0,sigma_a);
 }
 
 target+=log_sum_exp(to_vector(lp));
}

generated quantities{
   
     int<lower=1,upper=P*P> tau;
    // K simplex are a data type in Stan
    simplex[P*P] sp;
    
     vector[T] OD_pred;
     vector[T] DOD_pred;
     real MUGp;
     real MULp;
     real MUEp;
     
     MUGp=MUG;
     MULp=MUL;
     MUEp=MUE;
     sp = softmax(to_vector(lp));
     tau = categorical_rng(sp);
    for (n in 1:N){
      
      for (t in 1:T) {
        if(t<=L)
          OD_pred[t] =normal_rng(A_coef'*X[,t],sigma_o[1]);
        else{
          OD_pred[t] =normal_rng(A_coef'*X[,t],sigma_o[2]);
        }
        DOD_pred[t] = DA_coef'*X[,t];
      }
  }
}
