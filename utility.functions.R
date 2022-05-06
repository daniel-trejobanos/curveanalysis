library(rstan)
library(MASS)
model.code <-
'
functions {
  real normal_pdf(real obs, real mu, real sigma) {
    return exp(normal_log(obs,mu,sigma));
  }
  real normal_log_pdf(real obs,real mu, real sigma) {
    return normal_log(obs,mu,sigma);
  }
  real multi_normal_log_pdf(vector obs,vector mu, matrix SIGMA) {
    return multi_normal_log(obs,mu,SIGMA);
  }
  
}
model{}
'
expose_stan_functions(stanc(model_code = model.code))
#assuming sigma_o and sigma_d are known
mutual.information <- function(stan.fit,design,N,sigma_o,sigma_d) {
  M <- dim(design)[1]
  T <- dim(design)[2]

  acoef <- as.matrix(stan.fit,"A_coef")
  nsamples <- dim(acoef)[1]
  N1 <- nsamples/10
  N2 <- nsamples-N1
  IG <- rep(0,N1)
  SIGMA<-diag(rep(sigma_o^2,T))
  for(i in 1:N1){#sample from parameters 
    p <- 0
    X.q <- matrix(nrow=N,ncol=T)
    a.i <- matrix(acoef[i,],nrow=M,ncol=N)
    l.c <- 0
    for(n in 1:N){ #inner loop to evaluate per gene log
      X.q[n,] <- mvrnorm(1,mu=t(a.i[,n])%*%design,Sigma = SIGMA)  
      l.c <- l.c + multi_normal_log_pdf(X.q[n,],t(a.i[,n])%*%design,SIGMA)
    } 
    for(j in (N1+1):nsamples){
      log.c <-0
      for(n in 1:N){ #inner loop to evaluate per gene log
        a.j <- matrix(acoef[j,],nrow=M,ncol=N)
        log.c <- log.c + multi_normal_log_pdf(X.q[n,],t(a.j[,n])%*%design,SIGMA)
      } 
      p <- p + exp(log.c)
    }
    lp <- log(p/N2)
    IG[i] <-  l.c - lp
  }
  IG
}