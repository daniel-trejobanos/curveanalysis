add.iid.noise<-function(data,data.sd){
  N <- dim(data)[1]
  T <- dim(data)[2]
  data.set <- matrix(rnorm(N*T,mean = 0,sd = 0.02),N,T)
  data.out<- data+data.set  
}

subsample.data <- function(data,times,subsampling.rate){
  
  if(is.array(data)){
    N <- dim(data)[1]
    T <- dim(data)[2]
    if(subsampling.rate < T){
      subsamples.index <- seq(1,T,by = subsampling.rate)
      data.out=data[,subsamples.index]
    }
  }else{
    if(is.matrix(data)){
      N <- dim(data)[1]
      T <- dim(data)[2]
      if(subsampling.rate < T){
        subsamples.index <- seq(1,T,by = subsampling.rate)
        data.out[,subsamples.index]
      }
    }
  }
  list(data=data.out,times=times[subsamples.index])
}


