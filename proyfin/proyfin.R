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
datax<-list(x=c(1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1,
                  1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0,
                  0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1,
                  0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1),
            y= c(69, 57, 61, 60, 69, 74, 63, 68, 64, 53, 60, 58, 79, 56, 53, 74, 56, 76, 72,
                   56, 66, 52, 77, 70, 69, 76, 72, 53, 69, 59, 73, 77, 55, 77, 68, 62, 56, 68, 70, 60,
                   65, 55, 64, 75, 60, 67, 61, 69, 75, 68, 72, 71, 54, 52, 54, 50, 75, 59, 65, 60, 60,
                   57, 51, 51, 63, 57, 80, 52, 65, 72, 80, 73, 76, 79, 66, 51, 76, 75, 66, 75, 78, 70,
                   67, 51, 70, 71, 71, 74, 74, 60, 58, 55, 61, 65, 52, 68, 75, 52, 53, 70),
            z=c(1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0,
                   1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
                   1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1,1, 0, 1, 1, 0, 0, 1, 0, 0, 1),
            n=100)

inits<-function(){list(alpha=0, b.x=1, b.y=1)}
#-Selecting parameters to monitor-
parameters<-c("alpha","b.x","b.y","p")

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

