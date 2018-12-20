mu<-.4
sig<-.17
n<-10000
eta<-rnorm(n,mu,sig)
th<-exp(eta)/(1+exp(eta))
c(mean(th),sd(th))
