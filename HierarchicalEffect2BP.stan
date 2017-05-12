##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of coefficients
  int<lower=0> T; #number of timepoints
  int<lower=1> PL; #number of plates
  #int<lower=1> S;#number of strains
  int<lower=1,upper=PL> pl[N];#plate id
 # int<lower=1,upper=S> s[N];# strain id
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

parameters {
  vector<lower=0,upper=MAX-MIN>[M] A_coef;
  matrix[M,N] U[PL]; #well effects.
  real<lower=0> sigma_a;
  real<lower=0,upper=1> sigma_o[2,PL];
  
 real<lower=0> sigma_mu;
 real<lower=0> sigma_u;
 vector<lower=0>[3] lambda;
  real<lower=-MU[1],upper=0.01> VL;
  real<lower=-MU[2]> VG;
  real<lower=-MU[3]> VE;
}
transformed parameters{
    matrix[P,P] lp; 
    vector[M] DA_coef;
    row_vector[T] loss;
    row_vector[T] curve;
    row_vector[T] derivative;
    row_vector[T] logCurve;
    matrix[M,N] H_A_coef[PL];
      real MUL;
  real MUG;
  real MUE; 
    #vector<lower=0>[3] mu_strain;
   
    lp = rep_matrix(-log(P),P,P);
    for(i in 1:PL)
      H_A_coef[i]= rep_matrix(A_coef,N)+U[i];
    
    DA_coef=(1/MAXT)*(A_coef'*D)';
    curve=A_coef'*X;
    derivative=DA_coef'*X;
    loss = derivative ./ curve;
    logCurve=log(curve ./ curve[1]);
    MUL=MU[1]+VL;
    MUG=MU[2]+VG;
    MUE=MU[3]+VE;
   
   {
      vector[T+1] lpS_12;
      vector[T+1] lpS_21;
      vector[T+1] lpS_32;
      
      lpS_12[1]=0;
      lpS_21[1]=0;
      lpS_32[1]=0;
        for (t in 1:T){
            
            lpS_12[t+1]=lpS_12[t]+normal_lpdf(logCurve[t]|MUL,lambda[1]);
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
  sigma_o[1,]~ normal(0,0.1);
  sigma_o[2,]~ normal(0,0.1);
  sigma_a ~ normal(0,SIGMA_A);
  sigma_u~ cauchy(0,SIGMA_U);
   VL~normal(0,sigma_mu);
   VG~normal(0,sigma_mu);
   VE~normal(0,sigma_mu);
    A_coef[1]~cauchy(0,0.05);
    A_coef[2:M]  ~ normal(0,sigma_a);
 
 for(n in 1:N){
   U[pl[n],,n]~normal(0,sigma_u);
   OD[n,1:L]~normal(H_A_coef[pl[n],,n]'*X[,1:L],sigma_o[1,pl[n]]);
   OD[n,(L+1):T]~normal(H_A_coef[pl[n],,n]'*X[,(L+1):T],sigma_o[2,pl[n]]);
 }
 
   target+=log_sum_exp(to_vector(lp));
 
 
}

generated quantities{
     int s1;
     int Is2;
     int counter;
     real s2;
    
     real tauR;
     real PR;
     vector[T] OD_pred;
     vector[T] DOD_pred;
     vector[T] piece_pred;
     
     int<lower=1,upper=P*P> tau;
    simplex[P*P] sp;
    sp = softmax(to_vector(lp));
    tau = categorical_rng(sp);
    tauR=tau;
    PR=P;
     s1=modulus(tau,P);
     s2=ceil(tauR/PR);
     for(i in 1:P){
       if(i<=(s2+0.1))
          Is2=i;
     }
      for (t in 1:T) {
         OD_pred[t] = A_coef'*X[,t];
         DOD_pred[t] = DA_coef'*X[,t];
        if(t<prior[s1]){
          piece_pred[t]= OD_pred[1]*exp(MUL*MAXT*timeN[t]);
         
        }
        else{
          if(t<prior[Is2]){
           piece_pred[t]= OD_pred[prior[s1]]*exp(MUG*MAXT*(timeN[t]-timeN[prior[s1]]));
          
          }
          else
            piece_pred[t]= OD_pred[prior[Is2]]*exp(MUE*MAXT*(timeN[t]-timeN[prior[Is2]]));
        }
       
      }
  
}
