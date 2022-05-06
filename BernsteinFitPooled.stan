##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of coefficients
  int<lower=0> T; #number of timepoints
  #int<lower=1,upper=N> W #well id
   real MAXT;
  matrix<lower=0,upper=1>[M,T] X;
  matrix[M,M] D;
  real timeN[T];
  matrix[N,T] OD; #optical density data
   vector[3] MU;
  real MIN;
  real MAX;
  real SIGMA_MU;
  real SIGMA_A;
  real SIGMA_U;
   real<lower=0> LAMBDA;
    int P;
   int prior[P];
  int L;
   #real S;
}
transformed data{
  real<lower=0> MI;
  matrix[N,T] FPOD;
   MI=M;
  for(i in 1:N)
    for(j in 1:T)
      FPOD[i,j] = OD[i,j];
}
parameters {
  vector<lower=0,upper=MAX-MIN>[M] A_coef;
  matrix[M,N] U; #well effects.
  real<lower=0> sigma_a;
  real<lower=0,upper=1> sigma_o[2];
  
 real<lower=0> sigma_mu;
 real<lower=0> sigma_u;
 vector<lower=0>[3] lambda;
  vector[3] V;
}
transformed parameters{
    matrix[P,P] lp; 
    vector[M] DA_coef;
    row_vector[T] loss;
    row_vector[T] curve;
    row_vector[T] derivative;
    row_vector[T] logCurve;
    matrix[M,N] A_coef_hat;
      real MUL;
  real MUG;
  real MUE; 
    #vector<lower=0>[3] mu_strain;
   
    lp = rep_matrix(-log(P),P,P);
    
    A_coef_hat= rep_matrix(A_coef,N)+U;
    DA_coef=(1/MAXT)*(A_coef'*D)';
    curve=A_coef'*X;
    derivative=DA_coef'*X;
    loss = derivative ./ curve;
    logCurve=log(curve ./ curve[1]);
    #mu_strain= MU + V;
    MUL=MU[1]+V[1];
    MUG=MU[2]+V[2];
    MUE=MU[3]+V[3];
   
   {
      vector[T+1] lpS_12;
      vector[T+1] lpS_21;
      vector[T+1] lpS_32;
      
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
          lp[s1,s2]= lp[s1,s2] + lpS_32[T+1] + lpS_12[prior[s1]] +(lpS_21[prior[s2]]-lpS_21[prior[s1]])  - lpS_32[prior[s2]];
        }
      }
    
    
}

model { 
  lambda ~  cauchy(0,LAMBDA);
  sigma_mu ~ cauchy(0,SIGMA_MU);
  sigma_o~ normal(0,0.1);
  sigma_a ~ normal(0,SIGMA_A);
  sigma_u~ cauchy(0,SIGMA_U);
 # MUL~ normal(MU[1],sigma_mu);
  #MUG~ normal(MU[2],sigma_mu);
  #MUE~ normal(MU[3],sigma_mu);
   V~normal(0,sigma_mu);
    A_coef[1]~normal(0.07,0.01);
    A_coef[2:M]  ~ normal(0,sigma_a);
 
 for(n in 1:N){
   U[,n]~normal(0,sigma_u);
   FPOD[n,1:L]~normal(A_coef_hat[,n]'*X[,1:L],sigma_o[1]);
   FPOD[n,(L+1):T]~normal(A_coef_hat[,n]'*X[,(L+1):T],sigma_o[2]);
 }
 
   target+=log_sum_exp(to_vector(lp));
 
 
}

generated quantities{

     vector[T] OD_pred;
     vector[T] DOD_pred;
     int<lower=1,upper=P*P> tau;
    simplex[P*P] sp;
    sp = softmax(to_vector(lp));
    tau = categorical_rng(sp);
     
   #int<lower=1,upper=M> tau[N];
    #simplex[M] sp;
      for (t in 1:T) {
      
        OD_pred[t] = A_coef'*X[,t];
        DOD_pred[t] = DA_coef'*X[,t];
      }
  
}
