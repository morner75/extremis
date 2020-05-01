#### BS Basic Function ####

#' Penalty matrix
#'
#' @param n a scalar / length of random variable
#' @param d a scalar / order of difference
#' @return a matrix with n by n-d /  penalty matrix
#' @export
Kdmat<-function(n,d){
  Dd=cbind(0,diag(n-1))-cbind(diag(n-1),0)
  if (d==1) {return(t(Dd)%*%Dd)}
  else{ for(i in 2:d){
    Dd=(cbind(diag(n-d),0)-cbind(0,diag(n-d)))%*%Dd}}
  return(t(Dd)%*%Dd)}

#' Equally sapced knots
#'
#' @param x a numeric vector / a continuous random variable
#' @param knot.num a sclar / number of knots
#' @return a numeric vector / equally spaced knots
equiKnot<-function (x,knot.num){  ## equispace knot
  prob<-1:(knot.num)
  knots<-quantile(x,probs=prob/(length(prob)+1),type = 1)
  return(knots)}

#' Covariance matrix of coefficients of B-splines
#'
#' @param Bs a matrix / evaluations of x values from B-splines
#' @param sigma2 a scalar / variance of response y
#' @param tau2 a scalar / precision parameter for penalty
#' @param d scalar / degree of difference
#' @return a matrix  /  covariance matrix for beta
Sigmaf <- function(Bs,sigma2,tau2,d) solve((t(Bs)%*%Bs)/sigma2+Kdmat(ncol(Bs),d)/tau2)


#### Plotting  #######

#' Plotting Bayesian sample
#'
#' @param x a numeric vector / exploratory variable
#' @param y a numeric vector / response variable
#' @param true a numeric vector / true curve
#' @param trjs a matrix / trajectory matrix
#' @param beta a numeric vector / coefficients of B-splines
#' @param Bs a matrix / evaluations of x from B-splines
#' @param adjust a function / inverse of link funciton for plotting
#' @return posterior plot of meaan, 95% credible bands, along with y values (rug) and true curve (if applicable)
#' @export
# x=x ; y=y ; true=true ; beta=Out[,1:24]; Bs=Bs
PostPlot <- function(x,y,true,trjs=NULL,beta=NULL,Bs=diag(length(x)),adjust=NULL){
  if(is.null(trjs)) trjs <- beta%*%t(Bs)
  if(!is.null(adjust)) trjs <- adjust(trjs)
  post.mean <- colMeans(trjs)
  post.ci <- t(sapply(seq_along(x), function(i) quantile(trjs[,i],probs=c(0.025,0.975))))
  data <- data.frame(x=x,y=y,post.mean=post.mean,lower=post.ci[,1],higher=post.ci[,2])
  p <- ggplot(data=data)+
      {if(!is.null(y)) geom_point(aes(x=x,y=y),shape=16,alpha=0.5)}+
        geom_line(aes(x=x,post.mean))+
      {if(!is.null(true)) geom_line(aes(x=x,y=true),col="Red",linetype="dashed")}+
        geom_ribbon(aes(x=x,ymin=lower,ymax=higher),col="Gray",alpha=0.2)
  p
  }


