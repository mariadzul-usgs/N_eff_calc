

#mod = stan model output
#ind = index for parameter
#lag = lag for autocorrelation:
#warmup = number indicating how many iterations to discard
#Assumes 3 chains
require("rstan")
stan.get.rho<-function(mod,ind, lag, warmup){
  x<-as.vector(unlist(mod@sim$sample[[1]][ind]))[-c(1:warmup)]
  for(j in 2:3){x<-cbind(x,as.vector(unlist(mod@sim$sample[[j]][ind]))[-c(1:warmup)])}
  N<-dim(x)[1] #Number of iterations per chain
  M<-dim(x)[2] #Number of chains (assume 3)
  mu<-apply(x,2,mean)
  MU<-mean(mu)
  
  B<-(N/(M-1))*sum((mu-MU)^2)
  
  X<-numeric()
  S2<-numeric()
  for(j in 1:3){
    z<-numeric()
    s2<-numeric()
    for(i in 1:lag) {s2[i]<-(x[i,j]-mu[j])^2}
    for(i in (lag+1):N){
      z[i-lag]<-(x[i,j]-x[i-lag,j])^2
      s2[i]<-(x[i,j]-mu[j])^2}
    S2[j]<-(1/(N-1))*sum(s2)
    X[j]<-(1/(N-lag))*sum(z)}
  
  Vt<-(1/M)*sum(X)
  W<-(1/M)*sum(S2)
  
  var_hat<-((N-1)/N)*W+(1/N)*B
  rho_t<-1-(Vt/(2*var_hat))
  rho_t}



#This is the version in the Gelman et al. 2013 (Bayesian Data Analysis), 3rd Ed., pp: 286-287
#Stop calculating autocorrelation when the sum of two successive correlations is negative
#Disclaimer: Not sure why T has to be the first odd positive integer ?  Did not understand why/how T had to be odd.

N.eff.fun<-function(ind, mod, warmup){
  rho<-c(0,1)
  i<-0
  M<-3
  N<-length(as.vector(unlist(mod@sim$sample[[1]][1]))[-c(1:warmup)])
  while(sum(rho[(length(rho)-1):length(rho)])>0) {
    i<-i+1
    rho<-c(rho, stan.get.rho(mod,ind,i, warmup))}
  out<-(M*N)/(1+2*sum(rho[-c(1:2,(length(rho)-1):length(rho))]))
  out
}

#This is the version in the Stan manual
#Stop calculating autocorrelation when you get a negative autocorrelation:
N.eff.fun2<-function(ind, mod, warmup){
  rho<-1
  i<-0
  M<-3
  N<-length(as.vector(unlist(mod@sim$sample[[1]][1]))[-c(1:warmup)])
  while(rho[length(rho)]>0){
    i<-i+1
    rho<-c(rho, stan.get.rho(mod,ind,i, warmup))}
  if(length(rho)>2){out<-(M*N)/(1+2*sum(rho[-c(1,length(rho))]))}
  if(length(rho)==2){out<-(M*N)/sum(rho[-length(rho)])}
  out
}


#This is the output for the n_eff method described in BDA:
n.eff<-numeric()
for(j in 1:length(summary(MOD)$summary[,1])){n.eff[j]<-N.eff.fun(j, MOD, warmup)}
o<-match(names(MOD@sim$sample[[1]]),row.names(summary(MOD)$summary))
out<-cbind(round(n.eff),round(summary(MOD)$summary[o,9]))



#This is the output for the n_eff method described in Stan manual:
n.eff_Stan<-numeric()
for(j in 1:length(summary(MOD)$summary[,1])){n.eff_Stan[j]<-N.eff.fun2(j, MOD,warmup)}
o<-match(names(MOD@sim$sample[[1]]),row.names(summary(MOD)$summary))
out<-cbind(round(n.eff_Stan),round(summary(MOD)$summary[o,9]))

