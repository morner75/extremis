
#' Exploratory plot for discrete data
#'
#' @param data data
#' @return barchart
#' @export
Explot_disc <- function(data){
  dataw <- gather(data,key="Ex",value = "size",-x)
  ggplot(data=dataw)+
    geom_col(aes(x=x,y=size,fill=Ex),position = "dodge")
}



#' Posterior plot for discrete data
#'
#' @param supp a support of data
#' @param post posterior sample
#' @return line plot / histogram(density estimate)
#' @export
postPlot_disc <- function(supp,post,type=c("prob","variate")){
  ### for debugging
  # x=0:10 ; post=replicate(10,DPdisc(M=100,x=x),simplify=TRUE)
  ###
  type <- match.arg(type)
  dataw <- gather(data.frame(x=supp,post),key="traj",value = "values",-x)
  if(type=="prob"){
      ggplot(data=dataw)+
        geom_line(aes(x=x,y=values,col=traj),linetype="solid",size=0.1)+
        stat_summary(aes(x=x,y=values),fun.y="mean",geom="line",linetype="dashed",size=1)+
        theme(legend.position = "none")
    }else{
    ggplot(data=dataw,x=x,y=values)+
      geom_histogram(aes(x=values,y=..density..,fill=traj),alpha=0.5)+
      geom_density(aes(x=values,col=traj))+
      facet_wrap(vars(traj))
    }
  }


#' Posterior plot for DPM mixtures
#'
#' @param supp a support of data
#' @param post posterior sample
#' @return line plot / histogram(density estimate)
#' @export
postPlot_cont <- function(supp,post,type=c("mixture")){
  ### for debugging
  # supp=seq(-5,10,length.out = 200) ; post=Out1
  ###
  type <- match.arg(sample)
  dataw <- gather(data.frame(x=supp,t(post$trj(supp))),key="traj",value = "values",-x)
  if(type=="mixture"){
    ggplot(data=dataw)+
      geom_line(aes(x=x,y=values,col=traj),linetype="solid",size=0.1)+
      stat_summary(aes(x=x,y=values),fun.y="mean",geom="line",linetype="dashed",size=1)+
      theme(legend.position = "none")
  }
}
