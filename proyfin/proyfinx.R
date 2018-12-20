library(R2OpenBUGS)
library(R2jags)

wdir<-"/home/abraham/RA2018/proyfin"
setwd(wdir)

camp<-read.csv('multivsamp.csv',header=TRUE,sep=";")
campx<-camp[sample(nrow(camp), 1000), ]#camp[c(1:1000),]

n<-nrow(campx)

x<-campx$TO_NECESIDAD_FINAN_CAP_1M
y<-campx$NU_VINC_COGNODATA
z<-campx$contrata

#-Defining data-
data<-list("n"=n,"x"=x,"y"=y,"z"=z)

inits<-function(){list(alpha=0, b.x=1, b.y=1,zf=rep(0,n))}
#-Selecting parameters to monitor-
parameters<-c("alpha","b.x","b.y","zf")

#-Running code-
#OpenBUGS
ex1.sim<-bugs(data,inits,parameters,model.file="logreg.txt",
              n.iter=100000,n.chains=1,n.burnin=10000)
#OpenBUGS
out<-ex1.sim$sims.list

z<-out$b.x
par(mfrow=c(2,2))
plot(z,type="l")
plot(cumsum(z)/(1:length(z)),type="l")
hist(z,freq=FALSE)
acf(z)

out.sum<-ex1.sim$summary
print(out.sum)

traceplot(ex1.sim)
rex<-function(fun){return(fun)}

out<-rex(ex1.sim)$sims.list

out.sum<-rex(ex1.sim)$summary

#Predictions
adj<-function(out.sum,x,y,a){
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

adj(out.sum,y,z,"zf")

#DIC
met<-function(ex1.sim,out.sum){
  out.dic<-ex1.sim$DIC
  #out.dic<-ej4.sim$BUGSoutput$DIC
  print(out.dic)
  print(out.sum)
}

met(ex1.sim,out.sum)

r2<-function(y,out.sum,a){
  out.yf<-out.sum[grep(a,rownames(out.sum)),]  
  R2<-(cor(y,out.yf[1:n,1]))^2
  print(R2)
}

r2(y,out.sum,"zf")
