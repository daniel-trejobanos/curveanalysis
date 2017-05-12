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
   matrix[(M-1),T] XD;
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
   #real s;
}

parameters {
  vector<lower=0>[M] B;
  matrix[M,N] U[PL]; #well effects.
  real<lower=0> sigma_a;
  real<lower=0,upper=1> sigma_o[2,PL];
  
 real<lower=0> sigma_mu;
 real<lower=0> sigma_u;
 vector<lower=0>[3] lambda;
  real<lower=-MU[2]> VG;
  real<lower=-MU[3]> VE;
  vector<lower=0,upper=MAX-MIN>[M] A_coef;
  
    
}
transformed parameters{
    vector[P] lp; 
    vector[M-1] DA_coef;
    row_vector[T] loss;
    row_vector[T] curve;
    row_vector[T] derivative;
    row_vector[T] logCurve;
    matrix[M,N] H_A_coef[PL];
    
  real MUG;
  real MUE; 
    #vector<lower=0>[3] mu_strain;
    #A_coef[1]=B[1];
    #for(i in 2:M)
    #  A_coef[i]= A_coef[i-1]+s*B[i];
    lp = rep_vector(-log(P),P);
    for(i in 1:PL)
      H_A_coef[i]= rep_matrix(A_coef,N)+U[i];
    
    DA_coef=(1/MAXT)*(A_coef[2:M]-A_coef[1:(M-1)]);
    #DA_coef=(1/MAXT)*s*B[2:M];
    curve=A_coef'*X;
    derivative=DA_coef'*XD;
    loss = derivative ./ curve;
    logCurve=log(curve ./ curve[1]);
   
    MUG=MU[2]+VG;
    MUE=MU[3]+VE;
   
   {
      vector[T+1] lpS_12;
      vector[T+1] lpS_21;
      vector[T+1] lpS_32;
      
      lpS_12[1]=0;
      lpS_21[1]=0;
        for (t in 1:T){
            
            lpS_12[t+1]=lpS_12[t]+normal_lpdf(loss[t]|MUG,lambda[1]);
            lpS_21[t+1]=lpS_21[t]+normal_lpdf(loss[t]|MUE,lambda[2]);
        }
      
        for(s1 in 1:P)
        lp[s1]= lp[s1] + lpS_21[T+1] + lpS_12[prior[s1]] - lpS_21[prior[s1]];
        
      }
    
    
}

model { 
  lambda ~  cauchy(0,LAMBDA);
  sigma_mu ~ cauchy(0,SIGMA_MU);
  sigma_o[1,]~ normal(0,0.1);
  sigma_o[2,]~ normal(0,0.1);
  sigma_a ~ normal(0,SIGMA_A);
  sigma_u~ cauchy(0,SIGMA_U);
   VG~normal(0,sigma_mu);
   VE~normal(0,sigma_mu);
  
    A_coef[1]~cauchy(0,0.05);
     #B[2:M] ~chi_square(sigma_a/M);
    A_coef[2:M]  ~ normal(0,sigma_a);
 
 for(n in 1:N){
   U[pl[n],,n]~normal(0,sigma_u);
   OD[n,1:L]~normal(H_A_coef[pl[n],,n]'*X[,1:L],sigma_o[1,pl[n]]);
   OD[n,(L+1):T]~normal(H_A_coef[pl[n],,n]'*X[,(L+1):T],sigma_o[2,pl[n]]);
 }
 
   target+=log_sum_exp(to_vector(lp));
 
 
}

generated quantities{

     vector[T] OD_pred;
     vector[T] DOD_pred;
     int<lower=1,upper=P> tau;
    simplex[P] sp;
    sp = softmax(to_vector(lp));
    tau = categorical_rng(sp);
     
   #int<lower=1,upper=M> tau[N];
    #simplex[M] sp;
      for (t in 1:T) {
      
        OD_pred[t] = A_coef'*X[,t];
        DOD_pred[t] = DA_coef'*XD[,t];
      }
  
}
