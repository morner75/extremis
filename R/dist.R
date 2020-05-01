# 1. reprametrisation

# 2. link functions

# reciprocal
recip <- function(x) exp(1/x)
recipD1 <- function(x) -exp(1/x)/x^2

muRecip <- function(beta,Bs) recip(Bs%*%beta)
dRecip <- function(beta,Bs) recip.d1(Bs%*%beta)
vRecip <- function(beta,Bs) 1/recip(Bs%*%beta)^2
WRecip <- function(beta,Bs) diag(as.vector(dRecip(beta,Bs)^2/vRecip(beta,Bs)),nrow=nrow(Bs),ncol=nrow(Bs))
penSigmaRecip <- function(beta,Bs,tau2,d) solve(t(Bs)%*%WRecip(beta,Bs)%*%Bs+Kdmat(ncol(Bs),d)/tau2)
ycRecip <- function(beta,Bs,y) Bs%*%beta+diag(as.vector(1/dRecip(beta,Bs)))%*%(y-muRecip(beta,Bs))
proposalRecip <- function(beta,tau2) {   # IWLS proposal
  Sigmap <- penSigmaRecip(beta,Bs=Bs,tau2,d=1)
  meanp <- Sigmap%*%t(Bs)%*%WRecip(beta,Bs)%*%ycRecip(beta,Bs,y)
  betap <- as.vector(rmvnorm(n=1,mean= meanp, sigma = Sigmap ))
  list(betap=betap, meanp=meanp, Sigmap=Sigmap)
}


# inverse logit
invLogit <- function(x) exp(x)/(1+exp(x))
invLogitD1 <- function(x) exp(x)/(1+exp(x))^2

muLogit <- function(beta, Bs) invLogit(Bs%*%beta)
dLogit <- function(beta, Bs) invLogitD1(Bs%*%beta)
vLogit <- function(beta, Bs) invLogitD1(Bs%*%beta)
WLogit <- function(beta,Bs) diag(as.vector(dLogit(beta,Bs)^2/vLogit(beta,Bs)),nrow=nrow(Bs),ncol=nrow(Bs))
penSigmaLogit <- function(beta,Bs,tau2,d) solve(t(Bs)%*%WLogit(beta,Bs)%*%Bs+Kdmat(ncol(Bs),d)/tau2)
ycLogit <- function(beta,Bs,y) Bs%*%beta+diag(as.vector(1/dLogit(beta,Bs)))%*%(y-muLogit(beta,Bs))
proposalLogit <- function(beta,tau2) {   # IWLS proposal
  Sigmap <- penSigmaLogit(beta,Bs=Bs,tau2,d=1)
  meanp <- Sigmap%*%t(Bs)%*%WLogit(beta,Bs)%*%ycLogit(beta,Bs,y)
  betap <- as.vector(rmvnorm(n=1,mean= meanp, sigma = Sigmap ))
  list(betap=betap, meanp=meanp, Sigmap=Sigmap)
}


# 2. full conditional for y in different dist.
logJL <- function(beta, Bs, tau2, K, d, dist=c("Bernoulli","Exponential")){
  dist <- match.arg(dist)
  if(dist=="Bernoulli"){
    logL <- function(beta) sum(dbinom(y,1,prob=muLogit(beta, Bs),log=TRUE))
  }else if(dist=="Exponential"){
    logL <- function(beta) sum(dexp(y,rate=1/muRecip(beta, Bs),log=TRUE))
  }
  lbeta_tau2 <- function(beta,tau2,K,d) -(nrow(K)-d)/2*log(tau2)- t(beta)%*%K%*%beta/2/tau2
  logL(beta)+lbeta_tau2(beta,tau2,K,d)
}





