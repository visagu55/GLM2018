library(R2OpenBUGS)
library(R2jags)

wdir<-"/home/abraham/RA2018/ejemex/ex2018"
setwd(wdir)

#--- Funciones utiles ---
prob<-function(x){
  out<-min(length(x[x>0])/length(x),length(x[x<0])/length(x))
  out
}

#Leemos datos
bill<-read.table("http://allman.rhon.itam.mx/~lnieto/index_archivos/BillsMXc.csv",sep=",",header=TRUE)
bill
x<-bill$C#/1000
y<-bill$Y#/1000

library(tidyverse)
library(ggplot2)


#a)********************************************************************************************

datan<-gather(bill,val,I,x20:x500) %>%filter(.,I==1) %>%mutate(.,proporcion_falsos=Y/C)

explo <- ggplot(datan,aes(x=val,y=Y))+geom_boxplot()+xlab("val")+ylab("n")
exploanio <- ggplot(datan,aes(x=as.factor(Year),y=Y))+geom_boxplot()+xlab("anio")+ylab("n")
explogpo <- ggplot(datan,aes(x=val,y=Y))+geom_boxplot()+xlab("val")+ylab("n")+
  facet_wrap(~as.factor(Year), scale="free")



#plot(x,y,type = "p")
#plot(bill$Year,y)
#hist(bill$Year)
#cor(x,y)

#b)******************************************************************************************
n<-nrow(bill)
#-Defining data-
#data<-list(x1=sal$X1,x2=sal$X2,x3=sal$X3)
data<-list("n"=n,"y"=y,"e"=bill$C,"x50"=bill$x50,"x100"=bill$x100,"x200"=bill$x200,"x500"=bill$x500)
#-Defining inits-
inits<-function(){list(b=rep(0,5),yf=rep(0,n))}
parameters<-c("b","yf","prob")

#OpenBUGS
ex1<-bugs(data,inits,parameters,model.file="eja.txt",
               n.iter=10000,n.chains=2,n.burnin=1000)






res<-function(func,out){
  ex1.sim<-func
  #out<-ex1.sim$sims.list
  

 for(i in 1:5){   
  b<-out$b[,i]
  par(mfrow=c(1,1))
  plot(b,type="l",xlab = "beta")
  plot(cumsum(b)/(1:length(b)),type="l",xlab = "beta")
  hist(b,freq=FALSE, main = "beta")
  acf(b)
 }
  plot(out$b)
  
  #out.sum<-ex1.sim$summary
  #return(ex1.sim)
  
}

rex<-function(fun){return(fun)}

out<-rex(ex1)$sims.list

out.sum<-rex(ex1)$summary


res(ex1,out)


#Predictions
adj<-function(out.sum,y,a){
  out.yf<-out.sum[grep(a,rownames(out.sum)),]
  or<-order(x)
  ymin<-min(y,out.yf[,c(1,3,7)])
  ymax<-max(y,out.yf[,c(1,3,7)])
  par(mfrow=c(1,1))
  plot(x,y,ylim=c(ymin,ymax))
  lines(x[or],out.yf[or,1],lwd=2,col=2)
  lines(x[or],out.yf[or,3],lty=2,col=2)
  lines(x[or],out.yf[or,7],lty=2,col=2)
  
  plot(out.yf[1:n,1],y,type ="p")
  
}

adj(out.sum,y,"yf")

#DIC
met<-function(ex1.sim,out.sum){
  out.dic<-ex1.sim$DIC
  #out.dic<-ej4.sim$BUGSoutput$DIC
  print(out.dic)
  print(out.sum)
}

met(ex1,out.sum)

r2<-function(y,out.sum,a){
  out.yf<-out.sum[grep(a,rownames(out.sum)),]  
  R2<-(cor(y,out.yf[1:n,1]))^2
  print(R2)
}

r2(y,out.sum,"yf")

#d)*************************************************************************************************
out.sum[grep("prob",rownames(out.sum)),c(1,3,7)] 



#e)****************************************************************************************************

parameters<-c("b","yf","prob")
#OpenBUGS
ex2<-bugs(data,inits,parameters,model.file="eje.txt",
          n.iter=100000,n.chains=2,n.burnin=10000)
traceplot(ex2)


out2<-rex(ex2)$sims.list

out2.sum<-rex(ex2)$summary


res(ex2,out2)
#ajustes
adj(out2.sum,y,"yf")
#estimaciones y DIC
met(ex2,out2.sum)
#R2
r2(y,out2.sum,"yf")

out2.sum[grep("b",rownames(out2.sum)),c(1,3,7)] 

#f)***********************************************************************************************


#g)************************************************************************************************
out2.sum[grep("prob",rownames(out2.sum)),c(1,3,7)]  


#h)*************************************************************************************************
#El DIC de la logistica es 1297 vs. 1317 de la liga cloglog  y sus R2 se parecen 
#adj(out.sum,y,"yf")

out2.yf<-out2.sum[grep("yf",rownames(out2.sum)),]
yf<-out2.yf[1:n,1]
plot(yf,y,type ="p")
plot(yf,y,type ="p",col=bill$x20)
plot(yf,y,type ="p",col=bill$x50)
#i)************************************************************************************************
out2.sum[grep("prob",rownames(out2.sum)),c(1,3,7)] 
