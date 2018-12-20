library(R2OpenBUGS)
library(R2jags)
library(tidyverse)
library(forcats)
library(xtable)

wdir<-"/home/abraham/RA2018/proyfin/entregables"
setwd(wdir)

camp<-read.csv('segt-1.csv',header=TRUE,sep=",")

library(ggplot2)

tot<-camp %>% group_by(FH_REF)%>%summarise(total=sum(target))
ggplot(tot,aes(x=as.Date(FH_REF),y=total)) + geom_line(col="blue") +
  xlab("Fecha") +
  ylab("Creditos colocados")





ggplot(camp,aes(x=as.Date(FH_REF),y=target)) + geom_line(col="red") +facet_wrap(~TP_SEGMENTO_FINAL2)+
  xlab("Fecha") +
  ylab("Creditos colocados")

ggplot(camp, aes(x=fct_reorder(TP_SEGMENTO_FINAL2, target), y=target,fill=TP_SEGMENTO_FINAL2)) + 
  geom_boxplot() + coord_flip()


# Tabla resumen
camp %>% 
  group_by(TP_SEGMENTO_FINAL2) %>% 
  summarise(mediana = median(target),
            promedio = mean(target))  %>% 
  ungroup %>% 
  mutate(TP_SEGMENTO_FINAL2 = fct_reorder(TP_SEGMENTO_FINAL2, promedio)) %>% 
  arrange(desc(TP_SEGMENTO_FINAL2)) %>% 
  xtable(rownames = F, digits = 2)

c("préstamos", 
  median(camp$target),  
  mean(camp$target))


library(corrplot)
nums <- unlist(lapply(camp, is.numeric))  
corrplot(cor(camp[,nums]))

pairs(camp[,nums])

campx<-camp[camp$TP_SEGMENTO_FINAL2=='DIVORCIADO',]
campx.ts<-ts(campx,start=c(2016-08-01,1),end=c(2018-07-01,24),frequency=24)


n<-nrow(campx)
campx$Tiempo<-1:n

plot(campx.ts)
pairs(campx)
#cor(campx)

#x1<-campx$TO_PROM_TO_CARGOS_3M
x1<-scale(campx$TO_PROM_TO_CARGOS_3M)[1:n]
x2<-scale(campx$NU_VINC_COGNODATA)[1:n]
x3<-scale(campx$TO_NECESIDAD_FINAN_CAP_3M)[1:n]
x4<-scale(campx$IM_PROM_GASTOS_3M)[1:n]
x5<-scale(campx$IM_SUM_SDO_CORTE_1M)[1:n]
t<-campx$X
y<-campx$target

#-Defining data-

data<-list("n"=n,"y"=y,"x1"=x1,"x2"=x2,"x3"=x3,"x4"=x4,"x5"=x5)
#data<-list("n"=n,"y"=c(y[1:(n-6)],rep(NA,6)))

#inits<-function(){list(alpha=0,beta=rep(0,5),tau=1,yf1=rep(1,n))}
inits<-function(){list(alpha=rep(0,n),beta=matrix(0,nrow=5,ncol=n),tau=1,tau.a=1,tau.b=rep(1,5),yf1=rep(1,n))}

#inits<-function(){list(beta=rep(0,n))}
#inits<-function(){list(beta=rep(0,n),tau.b=1,yf1=rep(1,n))}
#inits<-function(){list(mu=rep(1,n),tau.b=1,yf1=rep(1,n))}
#-Selecting parameters to monitor-
#parameters<-c("beta","p")
parameters<-c("alpha","beta","tau","yf1")
#parameters<-c("tau.b","yf1","mu")

#-Running code-
#OpenBUGS
ex1.sim<-bugs(data,inits,parameters,model.file="logbin.txt",
              n.iter=10000,n.chains=1,n.burnin=1000)
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

#Tabla resumen
out.sum.t<-out.sum[grep("beta",rownames(out.sum)),c(1,3,7)]
out.sum.t<-cbind(out.sum.t,apply(out$beta,2,prob))
dimnames(out.sum.t)[[2]][4]<-"prob"
print(out.sum.t)

#DIC
out.dic<-ex1.sim$DIC
#out.dic<-ej8b.sim$BUGSoutput$DIC
print(out.dic)

#Predictions

out.yf<-out.sum[grep("yf1",rownames(out.sum)),]
y<-data$y
ymin<-min(y,out.yf[,c(1,3,7)])
ymax<-max(y,out.yf[,c(1,3,7)])


#x1 vs. y
x<-data$x1
par(mfrow=c(1,1))
plot(x,y,type="p",col="grey50",ylim=c(ymin,ymax))
points(x,out.yf[,1],col=2,pch=16,cex=0.5)
segments(x,out.yf[,3],x,out.yf[,7],col=2)
#x2 vs. y
x<-data$x2
par(mfrow=c(1,1))
plot(x,y,type="p",col="grey50",ylim=c(ymin,ymax))
points(x,out.yf[,1],col=2,pch=16,cex=0.5)

segments(x,out.yf[,3],x,out.yf[,7],col=2)
#x3 vs. y
x<-data$x3
par(mfrow=c(1,1))
plot(x,y,type="p",col="grey50",ylim=c(ymin,ymax))
points(x,out.yf[,1],col=2,pch=16,cex=0.5)
segments(x,out.yf[,3],x,out.yf[,7],col=2)

#x4 vs. y
x<-data$x4
par(mfrow=c(1,1))
plot(x,y,type="p",col="grey50",ylim=c(ymin,ymax))
points(x,out.yf[,1],col=2,pch=16,cex=0.5)
segments(x,out.yf[,3],x,out.yf[,7],col=2)

#x5 vs. y
x<-data$x5
par(mfrow=c(1,1))
plot(x,y,type="p",col="grey50",ylim=c(ymin,ymax))
points(x,out.yf[,1],col=2,pch=16,cex=0.5)
segments(x,out.yf[,3],x,out.yf[,7],col=2)

#t vs. y
x<-campx$Tiempo
par(mfrow=c(1,1))
plot(x,y,type="p",col="grey50",ylim=c(ymin,ymax))
points(x,out.yf[,1],col=2,pch=16,cex=0.5)
segments(x,out.yf[,3],x,out.yf[,7],col=2)
par(mfrow=c(1,1))
plot(x,y,type="l",col="grey50",ylim=c(ymin,ymax))
lines(x,out.yf[,1],col=2,cex=0.5)
lines(x,out.yf[,3],col=2,lty=2)
lines(x,out.yf[,7],col=2,lty=2)

#betas
out.beta<-out.sum[grep("beta",rownames(out.sum)),]
plot(out.beta[1:24,1],type="l")
plot(out.beta[105:208,1],type="l")
plot(out.beta[209:312,1],type="l")

#alpha
out.alpha<-out.sum[grep("alpha",rownames(out.sum)),]
plot(out.alpha[,1],type="l")

