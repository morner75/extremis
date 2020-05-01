#' Bayeian soomthing methods for binary responses
#'
#' @param x a numeric vector / covariate
#' @param y a numeric vector / response
#' @param no.iter a scalar / number of iteraions
#' @param burn.in a scalar / burn.in period
#' @param init a list of a scalar and a numeric vector / list of sigma2 and beta
#' @param Bs a matrix / evalutions of x from B-splines
#' @param knot a scalar / number of knots
#' @param d a scalar / order of difference for penalty
#' @return a list of MCMC sample from Gibbs sampling, posterior of trjectories and coefficients beta, sigma2 and tau2
#' @export
bpsBern <- function(y, x=NULL, no.iter=1500, burn.in=500, init, Bs, knot=20, d, sampler=c("Gibbs","IWLS")) {
  type <- match.arg(sampler)
  if(missing(d)) d=1
  if(missing(Bs)) Bs <- bSpline(x,knots=quantile(x,seq(0,1,length.out=knot+2)[2:(knot+1)]),order = 3,intercept = TRUE)
  if(missing(init)) {
    tau2 =1
    init <- initBer(tau2=tau2, y=y, Bs=Bs, d=d)
  }

  # set priors
  priors <- list(a0=0.001, b0= 0.001, a1=0.001, b1= 0.001)

  # initial values

  # penalty matrix
  K <- Kdmat(ncol(Bs),d)
  # iterations
  cat("Start time", format(Sys.time(), "%X"), "\n")  # start time
  if(sampler=="Gibbs"){
  tau2 <- tau2 ; beta=init$beta
  output <- samplerBernGibbs(N=no.iter,Bs=Bs,d=d,beta,tau2=tau2,priors=priors)
  } else{
  current <- list(betac = init$beta, meanc =init$beta, Sigmac = init$Sigma)
  # Define vector for storing simulated values (rows = iterations; columns = parameter)
  output <- matrix(0, no.iter,ncol(Bs)+1)
  accept.ratio <- numeric(length = no.iter)
  # Iterations
  res <- samplerBernIWLS(N=no.iter,Bs=Bs,d=d,current=current,tau2=tau2,K=K,priors = priors)
  output <- res$sample
  print(paste0("Acceptance ratio is ", mean(res$ratio)))
  }
  beta.post <-  output[(burn.in+1):no.iter,1:ncol(Bs)]
  tau2.post <- output[(burn.in+1):no.iter,ncol(Bs)+1]
  trjs <- beta.post%*%t(Bs)
  cat("Finish time", format(Sys.time(), "%X"), "\n") # Finish time
  list(trjs=trjs,beta=beta.post,tau2=tau2.post)
}


# Bernoulli response Gibbs sampler
samplerBernGibbs <- function(N,Bs,d,beta,tau2,priors){
  output <- matrix(0,nrow=N,ncol=ncol(Bs)+1)
  for(i in 1:N){
  li <- latent_update(y=y,mean=Bs%*%beta,sd=1)
  beta <- updateBeta(y=li,Bs,1,tau2,d=d)
  tau2 <- updateTau2(beta,priors$a1,priors$b1,d)
  output[i,] <- c(beta,tau2)
  if(i%%(N%/%5)==0) print(paste0(i,"th sampling done"))
  }
  return(output)
}

# Bernoulli response IWLS smapler
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


# initialisation
initBer <- function(tau2,y,Bs,d){
  K <- Kdmat(ncol(Bs),d)
  #init.beta <- function(y,Bs,K,tau2){
  lLikBer <- function(beta) sum(dbinom(y,1,prob=muLogit(beta, Bs),log=TRUE))
  scBer <- function(beta) sum((y-muLogit(beta, Bs))/(muLogit(beta, Bs)*(1-muLogit(beta, Bs))))
  lpenBer <- function(beta,tau2) lLikBer(beta)-t(beta)%*%K%*%beta/tau2/2
  scpenBer <- function(beta,tau2) scBer(beta)-K%*%beta/tau2
  out <- optim(rnorm(24),fn=function(beta,tau2) -lpenBer(beta,tau2), gr=function(beta,tau2) -scpenBer(beta,tau2), tau2=tau2,
               method = "BFGS",control=list(maxit=200),hessian = TRUE)
  list(beta=out$par, Sigma=solve(out$hessian))
}



# sampling function for Li
sampleLi <- function(y,beta) sapply(seq_along(y),function(i) (-1)^(y[i]+1)*abs(rnorm(1,mean=(Bs%*%beta)[i],1)))


