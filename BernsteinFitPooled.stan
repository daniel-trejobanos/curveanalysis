##stan model for bernstein polynomials
data {
  int<lower=0> N; #number of series
  int<lower=0> M; #number of coefficients
  int<lower=0> T; #number of timepoints
   real MAXT;
  matrix<lower=0,upper=1>[M,T] X;
  matrix[M,M] D;
  real timeN[T];
  matrix[N,T] OD; #optical density data
   real MU[3];
  real MIN;
  real MAX;
  real SIGMA_MU;
  real SIGMA_A;
   real<lower=0> LAMBDA;
    int P;
   int prior[P];
  
   #real S;
}
transformed data{
  real<lower=0> MI;
  matrix[N,T] FPOD;
   MI=M;
  for(i in 1:N)
    for(j in 1:T)
      #FPOD[i,j] = log(OD[i,j])-log(OD[i,1]);
      #FPOD[i,j] = log(OD[i,j])-MIN;
      FPOD[i,j] = OD[i,j];
}
parameters {
 #real<lower=0, upper=S> s;
  vector<lower=0,upper=MAX-MIN>[M] A_coef;
  real<lower=0> sigma_a;
  real<lower=0,upper=1> sigma_o;
#  real intercept[N];
  vector<lower=0>[3] rate;
 real<lower=0> sigma_mu;
 real<lower=0> lambda;
 # matrix<lower=0,upper=MAX>[M,N] jump;
#  real<lower=0,upper=MAX> k_0;
}
transformed parameters{
    matrix[P,P] lp[N]; 
    vector[M] DA_coef;
    vector[T] loss;
   # matrix<lower=0>[M,N] A_coef;
    for(n in 1:N){
    #  A_coef[,n] = cumulative_sum(s*jump[,n]);
      lp[n] = rep_matrix(-log(P),P,P);
    }
    DA_coef=(1/MAXT)*(A_coef'*D)';
    loss = ((DA_coef'*X)./(A_coef'*X))';
   
    for (s1 in 1:P) 
      for (s2 in 1:P) 
      for (t in 1:T){
        if(t<prior[s1]){
          lp[s1,s2] = lp[s1,s2] + normal_lpdf(loss[t]|rate[1],lambda);
          #exponential_lpdf((DA_coef[t-1,n]-rate[n,1])^2|LAMBDA);
        }else{
          if(t<prior[s2]){
            lp[s1,s2] = lp[s1,s2] + normal_lpdf(loss[t]|rate[2],lambda);
            #exponential_lpdf((DA_coef[t-1,n]-rate[n,2])^2|LAMBDA);
          }else{
            lp[s1,s2] = lp[s1,s2] + normal_lpdf(loss[t]|rate[3],lambda);
            #exponential_lpdf((DA_coef[t-1,n]-rate[n,3])^2|LAMBDA);
          }
        }
      } 
    
}

model { 
  lambda ~  normal(0,LAMBDA);
  sigma_mu ~ normal(0,SIGMA_MU);
  sigma_o~ normal(0,0.1);
  sigma_a ~ normal(0,SIGMA_A);
  rate[1]~ normal(MU[1],sigma_mu);
   rate[2]~ normal(MU[2],sigma_mu);
   rate[3]~ normal(MU[3],sigma_mu);
   #rate[,1]~ cauchy(0,SIGMA_MU);
   #rate[,2]~ cauchy(0,SIGMA_MU);
   #rate[,3]~ cauchy(0,SIGMA_MU);
 
    A_coef  ~ normal(0,sigma_a);
 
 
 for(n in 1:N){
   target+=log_sum_exp(to_vector(lp[n]))+normal_lpdf(FPOD[n,]|A_coef'*X,sigma_o);
    #target+=normal_lpdf(FPOD[n,]|A_coef'*X,sigma_o);
 }
 
}

generated quantities{
   
  # matrix[N,T] rate1_pred;
   # matrix[N,T] rate2_pred;
    # matrix[N,T] rate3_pred;
     vector[T] OD_pred;
     vector[T] DOD_pred;
     
     
   #int<lower=1,upper=M> tau[N];
    #simplex[M] sp;
   
      for (t in 1:T) {
      
     #    rate1_pred[n,t] = normal_rng(rate[n,1]*timeN[t],sigma_o[n]);
      #  rate2_pred[n,t] = normal_rng(rate[n,2]*timeN[t],sigma_o[n]);
      #  rate3_pred[n,t] =normal_rng(rate[n,3]*timeN[t],sigma_o[n]);
        OD_pred[t] = A_coef'*X[,t];
        DOD_pred[t] = DA_coef'*X[,t];
      }
  
}
