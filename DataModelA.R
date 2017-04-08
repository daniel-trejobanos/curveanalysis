library(deSolve)
.repression.nk<-function(x.kt,m.nk,theta.nk){
  1/((1+(x.kt/theta.nk)^m.nk))
}
.activation.nk<-function(x.kt,m.nk,theta.nk){
  (x.kt^m.nk)/(theta.nk^m.nk+x.kt^m.nk)
}

generate.data.model.A <- function(times){
  s<-rep(0.2,5)
  gamma<-c(0.9,0.9,0.7,1.5,1.5)
  m.nk=5
  theta.nk=1.5
  beta=2
  gene.expression<-function(t,y,params){
    dy1 <- s[1]-gamma[1]*y[1]+beta*.activation.nk(y[5],m.nk,theta.nk)
    dy2 <- s[2]-gamma[2]*y[2]+beta*.activation.nk(y[1],m.nk,theta.nk)
    dy3 <- s[3]-gamma[3]*y[3]+beta*.activation.nk(y[1],m.nk,theta.nk)
    dy4 <- s[4]-gamma[4]*y[4]+beta*.repression.nk(y[3],m.nk,theta.nk)+
      beta*.activation.nk(y[1],m.nk,theta.nk)
    dy5 <- s[5]-gamma[5]*y[5]+beta*.repression.nk(y[2],m.nk,theta.nk)+
      beta*.activation.nk(y[4],m.nk,theta.nk)
    list(c(dy1,dy2,dy3,dy4,dy5))       
  }
  
  yini<-c(1,0.5,1,0.5,0.5)
  
  out <- ode(func = gene.expression, parms = NULL, y = yini,
             times = times)
  matplot(out)
  N<-5
  T<-length(times)
  out<-out[,2:6]
  out<-t(out)
  
  node.1<-c(0,0,0,0,1) 
  node.2<-c(1,0,0,0,0)
  node.3<-c(1,0,0,0,0)
  node.4<-c(1,0,1,0,0)
  node.5<-c(0,1,0,1,0)
  H<-as.matrix(rbind(node.1,node.2,node.3,node.4,node.5))
  list(clean=out,adjacency=H)
}

generate.data.linear.model.A <- function(times){
  s<-rep(0.2,5)
  gamma<-c(0.9,0.9,0.7,1.5,1.5)
  m.nk=5
  beta=2
  gene.expression<-function(t,y,params){
    dy1 <- s[1]-gamma[1]*y[1]+beta*y[5]
    dy2 <- s[2]-gamma[2]*y[2]+beta*y[1]
    dy3 <- s[3]-gamma[3]*y[3]+beta*y[1]
    dy4 <- s[4]-gamma[4]*y[4]-beta*y[3]+
      beta*y[1]
    dy5 <- s[5]-gamma[5]*y[5]-beta*y[2]+
      beta*y[4]
    list(c(dy1,dy2,dy3,dy4,dy5))       
  }
  
  yini<-c(1,0.5,1,0.5,0.5)
  
  out <- ode(func = gene.expression, parms = NULL, y = yini,
             times = times)
  matplot(out)
  N<-5
  T<-length(times)
  out<-out[,2:6]
  out<-t(out)
  
  node.1<-c(0,0,0,0,1) 
  node.2<-c(1,0,0,0,0)
  node.3<-c(1,0,0,0,0)
  node.4<-c(1,0,1,0,0)
  node.5<-c(0,1,0,1,0)
  H<-as.matrix(rbind(node.1,node.2,node.3,node.4,node.5))
  list(clean=out,adjacency=H)
}