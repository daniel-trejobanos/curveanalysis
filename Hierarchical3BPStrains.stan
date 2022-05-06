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
  matrix<lower=0>[M-1,T] XD;
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
  int LAGT;
  
  #real S;
}
transformed data{
  matrix[M,M] identity; //Identity for scaling matrix

  identity = diag_matrix(rep_vector(1.0,M)); 
}
parameters {
  matrix<lower=0>[M,PL] A_coef;
  #matrix[M,N] U[PL]; #well effects.
  #real<lower=0> sigma_a;
  matrix<lower=0,upper=1>[2,PL] sigma_o;
  
  #matrix<lower=0>[3,PL] sigma_mu;
  #real<lower=0> sigma_u;
  vector<lower=0>[4] lambda;
  real<lower=0,upper=0.001> LAG[PL];
  real<lower=-MU[1],upper=1> VL[PL];
  real<lower=-MU[2],upper=1> VG[PL];
  real<lower=-MU[3],upper=1> VE[PL];
  cov_matrix[M] SIGMA[PL];
  matrix<lower=0,upper=0.2>[M,PL] alpha;
#  matrix<lower=0,upper=1>[2,PL] sigma_o_raw;
}
transformed parameters{
  matrix[P,P] lp[PL]; 
  #matrix[M,PL] DA_coef;
  matrix[M-1,PL] DA_coef;
 matrix[PL,T] loss;
  matrix[PL,T] curve;
  matrix[PL,T] derivative;
 # matrix[PL,T] logCurve;
 # matrix[M,N] H_A_coef[PL];
  real MUL[PL];
  real MUG[PL];
  real MUE[PL]; 
  real MULAG[PL];
  
  
  matrix<lower=0>[M,PL] betaO;
  #vector<lower=0>[3] mu_strain;
   for(i in 1:PL)
   lp[i] = rep_matrix(-log(P),P,P);
   
   #sigma_o[1,]=0.001*sigma_o_raw[1,];
   #sigma_o[2,]=0.001*sigma_o_raw[2,];
                    
  for(i in 1:PL){
    betaO[1,i]=alpha[1,i];
     for(k in 2:M)
        betaO[k,i]=betaO[k-1,i]+alpha[k,i];
    
  
   DA_coef[,i]=(1/MAXT)*(A_coef[2:M,i]-A_coef[1:(M-1),i]);
                    curve[i,]=A_coef[,i]'*X;
                    derivative[i,]=DA_coef[,i]'*XD;
                    loss[i,] = derivative[i,] ./ curve[i,];
 #   logCurve[i,]=log(curve[i,] ./ curve[i,1]);
   # H_A_coef[i]= rep_matrix(A_coef[,i],N)+U[i];
  
 
                    MUL[i]=MU[1]+VL[i];
                    MUG[i]=MU[2]+VG[i];
                    MUE[i]=MU[3]+VE[i];
                    MULAG[i]=0+LAG[i];
                    {
                      vector[T+1] lpS_12;
                      vector[T+1] lpS_21;
                      vector[T+1] lpS_32;
                      vector[T+1] lpS_lag;
                      
                      lpS_12[1]=0;
                      lpS_21[1]=0;
                      lpS_32[1]=0;
                       lpS_lag[1]=0;
                      for (t in 1:T){
                        lpS_lag[t+1]=lpS_lag[t]+normal_lpdf(loss[i,t]|MULAG[i],lambda[1]);
                        lpS_12[t+1]=lpS_12[t]+normal_lpdf(loss[i,t]|MUL[i],lambda[2]);
                        lpS_21[t+1]=lpS_21[t]+normal_lpdf(loss[i,t]|MUG[i],lambda[3]);
                        lpS_32[t+1]=lpS_32[t]+normal_lpdf(loss[i,t]|MUE[i],lambda[4]);
                      }
                      
                      for(s1 in 1:P-1)
                        for(s2 in (s1+1):P){
                          lp[i,s1,s2]= lp[i,s1,s2] + lpS_32[T+1] + lpS_lag[LAGT]+ (lpS_12[prior[s1]]-lpS_12[LAGT]) +(lpS_21[prior[s2]]-lpS_21[prior[s1]])  - lpS_32[prior[s2]];
                        }
                    }
    
  }
                    
                    
}

model { 
  
  lambda ~  cauchy(0,LAMBDA);
  
  sigma_o[1,]~ normal(0,0.001);
  sigma_o[2,]~ normal(0,0.01);
  #sigma_a ~ normal(0,SIGMA_A);
  #sigma_u~ cauchy(0,SIGMA_U);
 
  #for(i in 1:PL){
   #A_coef[1,i]~cauchy(0,0.05);
   #A_coef[2:M,i]  ~ normal(0,sigma_a);
    #target+=log_sum_exp(to_vector(lp[i]));
  #}
  
 
  for(n in 1:N){
    alpha[1,pl[n]]~normal(0,0.1);
    alpha[2:M,pl[n]]~normal(0,1);
    
    #A_coef[1,pl[n]]~normal(0,0.1);
    #A_coef[2:M,pl[n]]  ~ normal(0,sigma_a);
   # U[pl[n],,n]~normal(0,sigma_u);
   #sigma_mu[,pl[n]] ~ cauchy(0,SIGMA_MU);
   # VL[pl[n]]~normal(0,sigma_mu[1,pl[n]]);
    #VG[pl[n]]~normal(0,sigma_mu[2,pl[n]]);
    #VE[pl[n]]~normal(0,sigma_mu[3,pl[n]]);
    VL[pl[n]]~normal(0,0.3);
    VG[pl[n]]~normal(0,0.1);
    VE[pl[n]]~normal(0,0.01);
    
   SIGMA[pl[n]] ~inv_wishart(M+2,0.01*identity);
   A_coef[,pl[n]] ~ multi_normal(betaO[,pl[n]],SIGMA[pl[n]]);
    OD[n,1:L]~normal(A_coef[,pl[n]]'*X[,1:L],sigma_o[1,pl[n]]);
   OD[n,(L+1):T]~normal(A_coef[,pl[n]]'*X[,(L+1):T],sigma_o[2,pl[n]]);
   target+=log_sum_exp(to_vector(lp[pl[n]]));
  }
  
 
  
  
}

generated quantities{
     #int s1;
     #int Is2;
     #int counter;
     #real s2;
     #real tauR;
     #real PR;
     vector[T] OD_pred;
     vector[T] DOD_pred;
     matrix[N,T] log_lik;
    # vector[T] piece_pred;
     
     int<lower=1,upper=P*P> tau;
    simplex[P*P] sp;
    sp = softmax(to_vector(lp[1]));
    tau = categorical_rng(sp);
    #tauR=tau;
    #PR=P;
    # s1=modulus(tau,P);
    # s2=ceil(tauR/PR);
     #for(i in 1:P){
    #   if(i<=(s2+0.1))
    #      Is2=i;
    # }
      for (t in 1:T) {
         OD_pred[t] = A_coef[,1]'*X[,t];
         DOD_pred[t] = DA_coef[,1]'*XD[,t];
         for(n in 1:N){
            if(t<=L)
              log_lik[n,t]=normal_lpdf(OD[n,t]|A_coef[,pl[n]]'*X[,t],sigma_o[1,pl[n]]);
            else      
              log_lik[n,t]=normal_lpdf(OD[n,t]|A_coef[,pl[n]]'*X[,t],sigma_o[2,pl[n]]);  
         }
         
     #    if(t<LAGT){
      #     piece_pred[t]= OD_pred[1]*exp(MULAG[1]*MAXT*timeN[t]);
       #  }
        #else{
         # if(t<prior[s1]){
         # piece_pred[t]= OD_pred[LAGT]*exp(MUL[1]*MAXT*(timeN[t]-timeN[LAGT]));
         
        #}
        #else{
         # if(t<prior[Is2]){
          # piece_pred[t]= OD_pred[prior[s1]]*exp(MUG[1]*MAXT*(timeN[t]-timeN[prior[s1]]));
          
          #}
          #else
           # piece_pred[t]= OD_pred[prior[Is2]]*exp(MUE[1]*MAXT*(timeN[t]-timeN[prior[Is2]]));
        #}
        #}
          
       
      }
  
}

