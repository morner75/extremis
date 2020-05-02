#### BS Basic Function ####
Kdmat<-function(n,d){
  Dd=cbind(0,diag(n-1))-cbind(diag(n-1),0)
  if (d==1) {return(t(Dd)%*%Dd)}
  else{ for(i in 2:d){
    Dd=(cbind(diag(n-d),0)-cbind(0,diag(n-d)))%*%Dd}}
  return(t(Dd)%*%Dd)}

# Equally sapced knots
equiKnot<-function (x,knot.num){  ## equispace knot
  prob<-1:(knot.num)
  knots<-quantile(x,probs=prob/(length(prob)+1),type = 1)
  return(knots)}

# Covariance matrix of coefficients of B-splines
Sigmaf <- function(Bs,sigma2,tau2,d) solve((t(Bs)%*%Bs)/sigma2+Kdmat(ncol(Bs),d)/tau2)




plot.chit <- function(res,true=NULL){
  trjs <- pnorm(res$trjs)
  post.mean <- colMeans(trjs)
  post.ci <- t(sapply(seq_along(res$time), function(i) quantile(trjs[,i],probs=c(0.025,0.975))))
  data <- data.frame(x=res$time,y=res$y,post.mean=post.mean,lower=post.ci[,1],higher=post.ci[,2])
  p <- ggplot(data=data)+
    geom_point(aes(x=x,y=y),shape=16,alpha=0.5)+
    geom_line(aes(x=x,post.mean))+
    {if(!is.null(true)) geom_line(aes(x=x,y=true),col="Red",linetype="dashed")}+
    geom_ribbon(aes(x=x,ymin=lower,ymax=higher),col="Gray",alpha=0.2)
  p
}

