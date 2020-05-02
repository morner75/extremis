

# Time-varying extremal coefficient chi function using Bayeian soomthing methods for binary responses

# S3 methods
chit <- function(It, time, no.iter, burn.in=500, sampler=c("Gibbs","IWLS"),...)
  UseMethod("chit")


# Sampling chit
chit.default <- function(y, time, no.iter=1500, burn.in=500, sampler=c("Gibbs","IWLS"), init, knot=20, d=1) {
  type <- match.arg(sampler)
  Bs <- bSpline(time,knots=quantile(time,seq(0,1,length.out=knot+2)[2:(knot+1)]), order = 3,intercept = TRUE)

  # initial values
  if(missing(init)){
    tau2=1
    init <- initBer(tau2=tau2, y=y, Bs = Bs, d=d)
  }

  # set priors
  priors <- list(a0=0.001, b0= 0.001, a1=0.001, b1= 0.001)
  
  # penalty matrix
  K <- Kdmat(ncol(Bs),d)
  
  # start time
  cat("Start time", format(Sys.time(), "%X"), "\n")  

  # call sampler
  if(sampler=="Gibbs"){  # Gibbs
    
    tau2 <- tau2 ; beta=init$beta
    out <- samplerBernGibbs(N=no.iter,y=y,Bs=Bs,d=d,beta,tau2=tau2,priors=priors)
    
  } else{ # IWLS
    
    current <- list(betac = init$beta, meanc =init$beta, Sigmac = init$Sigma)
    # Define vector for storing simulated values (rows = iterations; columns = parameter)
    out <- matrix(0, no.iter,ncol(Bs)+1)
    accept.ratio <- numeric(length = no.iter)
    # Iterations
    res <- samplerBernIWLS(N=no.iter,Bs=Bs,d=d,current=current,tau2=tau2,K=K,priors = priors)
    out <- res$sample
    print(paste0("Acceptance ratio is ", mean(res$ratio)))
    
  }
  beta.post <-  out[(burn.in+1):no.iter,1:ncol(Bs)]
  tau2.post <- out[(burn.in+1):no.iter,ncol(Bs)+1]
  trjs <- beta.post%*%t(Bs)
  cat("Finish time", format(Sys.time(), "%X"), "\n") # Finish time
  outputs <- list(y=y,time=time,trjs=trjs,beta=beta.post,tau2=tau2.post)
  class(outputs) <- "chit"
  return(outputs)
  }


# Gibbs sampler
samplerBernGibbs <- function(N,y,Bs,d,beta,tau2,priors){
  output <- matrix(0,nrow=N,ncol=ncol(Bs)+1)
  for(i in 1:N){
  li <- updateLatent(y=y,beta=beta,Bs=Bs)
  beta <- updateBeta(y=li,Bs,1,tau2,d=d)
  tau2 <- updateTau2(beta,priors$a1,priors$b1,d)
  output[i,] <- c(beta,tau2)
  if(i%%(N%/%5)==0) print(paste0(i,"th sampling done"))
  }
  return(output)
}

# IWLS smapler
samplerBernIWLS <- function(N,Bs,d,current,tau2,K,priors){
  output <- matrix(0, N,ncol(Bs)+1)
  alpha <- accept.ratio <- numeric(length = N)
  for (i in 1:N) {
    proposal <- proposalLogit(current$betac,tau2)
    beta.update <- updateBetaBernIWLS(proposal,current,tau2,K,d)
    current <- beta.update$current
    accept.ratio[i] <- beta.update$status
    tau2 <- updateTau2(beta=current$betac,a1=priors$a1,b1=priors$b1,d=1)
    output[i, ] <- c(current$betac,tau2)
    if(i%%(N%/%5)==0) print(paste0(i,"th sampling done"))
  }
  m.ratio <- mean(accept.ratio)
  return(list(sample=output,ratio=m.ratio))
}

# sampling function for beta with IWLS proposals
updateBetaBernIWLS <- function(proposal,current,tau2,K, d){
  current.likelihood <- logJL(current$betac, Bs, tau2, K, d, dist = "Bernoulli")
  proposal.likelihood <- logJL(proposal$betap, Bs, tau2,K, d, dist = "Bernoulli")
  proposal.eval <- dmvnorm(proposal$betap,mean=current$meanc,sigma=current$Sigmac,log = TRUE)
  current.eval <- dmvnorm(current$betac, mean=proposal$meanp, sigma=proposal$Sigmap, log = TRUE)
  alpha <- proposal.likelihood-current.likelihood - proposal.eval  + current.eval
  if(exp(alpha) > runif(1)){
    list(current=list(betac=proposal$betap,meanc=proposal$meanp, Sigmac=proposal$Sigmap),alpha=exp(alpha),status=1)
  } else {
    list(current=current,alpha=exp(alpha),status=0)}
}

# update latent y for Gibbs
updateLatent <- function(y,beta,Bs){
  mean <- Bs%*%beta
  res = numeric(length(y))
  for(i in seq_along(res)){
    if(y[i]==1){
      while(res[i]<=0) res[i] = rnorm(1,mean[i],1)
    } else{
      while(res[i]>=0) res[i] = rnorm(1,mean[i],1)
    }
  }
  return(res)
}

# update beta for Gibbs
updateBeta <- function(y,Bs,sigma2,tau2,d) as.vector(rmvnorm(1, mean=Sigmaf(Bs,sigma2,tau2,d)%*%t(Bs)%*%y/sigma2, 
                                                             sigma= Sigmaf(Bs,sigma2,tau2,d)))
# update tau2 for Gibbs and IWLS
updateTau2 <- function(beta,a1,b1,d) 1/rgamma(1,a1+(length(beta)-d)/2,b1+t(beta)%*%Kdmat(length(beta),d)%*%beta/2)


# initialisation
initBer <- function(tau2,y,Bs,d){
  neglpenBer <- function(beta, y, tau2, Bs) -sum(dbinom(y,1,prob=muLogit(beta, Bs),log=TRUE)) + t(beta)%*%K%*%beta/tau2/2
  negscpenBer <- function(beta, y, tau2, Bs) -sum((y-muLogit(beta, Bs))/(muLogit(beta, Bs)*(1-muLogit(beta, Bs)))) + K%*%beta/tau2
  K <- Kdmat(ncol(Bs),d)
  out <- optim(rnorm(ncol(Bs)),fn=neglpenBer, gr=negscpenBer, y=y, tau2=tau2, Bs=Bs, method = "BFGS",control=list(maxit=200),hessian = TRUE)
  list(beta=out$par, Sigma=solve(out$hessian))
}

