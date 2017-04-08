bernstein.basis <- function(m,t){
  k=1:m
  choose(m-1,k-1)*(t^(k-1))*((1-t)^(m-k))
}
  
bernstein.matrix <- function(m,time){
  sapply(X = time,FUN=function(x)bernstein.basis(m,x))
}

bernstein.basis.derivative<- function(m,t){
  (m-1)*BernsteinBasis(m-1,t)
}

bernstein.basis.derivative.matrix<- function(m,time){
  sapply(X = time,FUN=function(x)bernstein.basis.derivative(m,x))
}


normalize.interval <- function(x,a=x[1],b=x[length(x)]){
  (x-a)/(a-b)
}

derivative.coefficients <- function(coefficients,order=1){
  M<-length(coefficients)
  if(order==1){
   (M-1)*(coefficients[2:length(coefficients)]-coefficients[1:(length(coefficients)-1)])
  }
  else{
    coef<-(coefficients[2:length(coefficients)]-coefficients[1:(length(coefficients)-1)])
    (M-1)*(derivative.coefficients(coef,order=order-1))
  }
}



get.roots<- function(coefficients){
  L <- length(coefficients)
  companion<- as.matrix(cbind(diag(L-2),-coefficients[2:(L-1)]))
  companion <- rbind(c(rep(0,L-2),-coefficients[1]),companion)
  eigen(companion)$values
}

coef.bernstein.to.power <- function(x){
  m <- length(x)
  k <- 1:m
  result<-x*choose(m-1,k-1)
  #the companion matrix requires the first coefficient to be 1
  if(result[m]!=0)
    result <- result/result[m]
  else{
    result <- result[1:m-1]
    result <- result/result[m-1]
  }
  result
}


get.maximum.values <- function(x){
  L<-length(x)
  m <- L
  dx <- derivative.coefficients(x) 
  cdx<-coef.bernstein.to.power(dx)
  tmp.t<-get.roots(cdx)
  tmp.t<-tmp.t/(1+tmp.t)
  tmp.t <- tmp.t[abs(Im(tmp.t)) < 1e-6]
  tmp.t<- Re(tmp.t)
  time<-tmp.t
  X <-BernsteinMatrix(m,time)
  tmp.v<-x%*%X
  
  value<-t(tmp.v)
  
   result<- cbind(t(t(as.array(time))),value,deparse.level = 0)
   result<-result[result[,2]==max(value),]
   result
}


bernstein.derivative.operator<- function(m){
  n.row = m-1
  n.col = m
  D <- matrix(nrow = n.row,ncol = n.col)
  for(i in 1:n.row){
    for(j in 1:n.col){
      if(i==j){
        D[i,j] <- -1
      }else{
        if(j==(i+1)){
          D[i,j] <- 1
        }else{
          D[i,j] <- 0
        }
        
      }
      
    }
  }
  L <- matrix(ncol=n.col,nrow=n.row)
  for(i in 1:n.row){
    for(j in 1:n.col){
      L[i,j]=0
      if(i==j){
        L[i,j]=((n.col-1)-(i-1))/(n.col-1)
      }
      if(j==(i+1)){
        L[i,j]=(j-1)/(n.col-1)  
      }
    }
  }
  
  (m-1)*t(D)%*%L
}