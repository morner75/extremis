bevfd.kernel<- function(Y, Xt, X, tau = 0.95, smoothing  = "lscv", 
                        grid = seq(0.01, 0.99, length = 2^8),
                        raw = TRUE, target = FALSE)
  UseMethod("bevfd.kernel")


bevfd.kernel.default<- function(Y, Xt, X, tau = 0.95, smoothing  = "lscv",
                          grid = seq(0.01, 0.99, length = 2^8), 
                          raw = TRUE, target = "pdf"){
    ## Basic input validation
    if (tau >= 1)
        stop('tau must be smaller than 1')
    if (dim(Y)[2] > 2)
        stop('this function only applies to bivariate extremes')
    if (target != "pdf" & target != "cdf" & target != "Pickands")
    stop('target must be equal to "pdf", "cdf", or "Pickands"')
     n <- dim(Y)[1]
    ## Convert data to common margins, if raw == TRUE
    if (raw == TRUE) {
        FX <- ecdf(Y[, 1])
        FY <- ecdf(Y[, 2])
        x <- -1 / log(n / (n + 1) * FX(Y[, 1]))
        y <- -1 / log(n / (n + 1) * FY(Y[, 2]))
    }
    if (raw == FALSE) {
        x <- Y[, 1]
        y <- Y[, 2]
    }
    ## Compute angles and threshold radius
    w <- x / (x + y)
    u <- quantile(x + y, tau)     
    w <- w[which(x + y > u)]
    
    if ( is.list( smoothing ) == TRUE ){
      
      b <- smoothing$b; nu <- smoothing$nu; c <- smoothing$c
      
     }else if( smoothing == "lscv" ){
       
      lscv<-function( nobreaks, smoothing ){
       b <- smoothing[1]; nu <- smoothing[2]; c <- smoothing[3]
        if(b <= 0)
          stop('b must be positive')
        if(nu <= 0)
          stop('nu must be positive')
        if (c < 0)
          stop('c must be positive')
       
        index<- 1 : length(w)
        s<- floor( length(w) / nobreaks )
        scv<-numeric()
        for(nk in 1 : nobreaks){
          ind<- index[ ( s * ( nk - 1) + 1) : ( s* nk) ]
          nm <- fda.usc::metric.lp(X, Xt)[1, ][ 1 : length( w[-ind] )]
          p <- 2 * dnorm(nm/b) / b * (nm > 0) / 
            sum(2 * dnorm(nm / b) / b * (nm > 0))
          h1<-NULL
          for( q in 1:length(w)){
            hestimator <- function(W){
              h<- sum( p * dbeta( W, (nu * w[-ind][q] * ( (.5) / 
            ( sum(p * w[-ind][q])))) + c,
            ( nu * (1 - (w[-ind][q] * (.5 / ( sum(p * w[-ind][q] )))))) + c) )
              h^2
            }
            W = seq(0.01, 0.99, length.out = 2^8)
            temp = tryCatch(sintegral( W, hestimator )$value, 
                            error = function(e) e)
            if(!inherits(temp, "error")){
              h1[q] = temp
            }else{
              h1[q] = 1e09
            }
            
          }
          hestimator <- function(W) sum( p * dbeta( W, (nu * w[-ind] * ( (.5) / 
                 ( sum(p * w[-ind])))) + c,
                 ( nu * (1 - (w[-ind] * (.5 / ( sum(p * w[-ind] )))))) + c) )
          
          h2 <- sapply( w[ind], hestimator)
          scv[nk]<- mean(h1) - 2*mean(h2)
        }
        return( sum( scv ) )
      }
      optimls<-optim( lscv, nobreaks = 10, par = list(1,1,14) )$par
      b <- optimls[1]; nu <- optimls[2]; c <- optimls[3]
      
   }
   nm <- fda.usc::metric.lp(X, Xt)[1, ][ 1 : length(w) ]
    p <- 2 * dnorm(nm/b) / b * (nm > 0) / 
        sum(2 * dnorm(nm / b) / b * (nm > 0))
      
    hestimator <- function(W) sum(p * dbeta( W, (nu * w * ( (.5) / 
    ( sum(p * w)))) + c, ( nu * (1 - (w * (.5 / (sum(p * w) ))))) + c) )
      
      h <- sapply(grid, hestimator)
   
  
  
    if (target == "pdf") {
      hestimator <-function(W) sum(p * dbeta(W, ( nu * w * ( (.5) / 
        ( sum( p * w )))) + c, ( nu * (1 - (w * (.5 / (sum( p * w) ))))) + c))
        h <- sapply( grid, hestimator)
        }
    
    if (target == "cdf") {
      hestimator <-function(W) sum( p * pbeta(W, (nu * w * ((.5) / 
     ( sum( p * w )))) + c, ( nu * (1 - (w * (.5 / (sum( p * w) ))))) + c))
      h <- sapply( grid, hestimator)
    }
     
    if (target == "Pickands") {
      k=10
    hestimator <- function(W) {
      ptilde <- numeric(k)
      integrand <- function(W) pbeta(W, (nu * w[i] * ((0.5) / 
     (sum(p * w[i])))) + c,(nu * (1 - (w[i] * (0.5 / (sum(p * w[i])))))) + c)
      for(i in 1:k) 
      ptilde[i] <- p[i] * as.numeric(integrate(integrand, 0, W)[1]) 
      1 - W + 2 * sum(ptilde)
      }
         hestimator <- Vectorize(hestimator)
         h <- sapply(grid, hestimator)
     }
    
    ## Organize and return outputs    
    outputs <- list(h = h, w = w, grid = grid,
               smoothing = list(b=b, nu=nu, c=c), target = target)
    class(outputs) <- "bevfd.kernel"
    return(outputs)
}


plot.bevfd.kernel<-  function(fit, xlab = "w", ...) {
  grid <- fit$grid
  trajhat <- fit$h
  target <- fit$target
  
  if (target == "pdf" | target == "cdf") {
    frame <- data.frame(grid = grid, trajhat = trajhat)
    gg <- ggplot(data = frame, aes(x = grid, y = trajhat)) +
      geom_line(aes(x = grid, y = trajhat), size = 1.1, col = "blue") +
      theme_classic() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 10))
    if (target == "pdf")
      gg <- gg + labs(y = "Cross Section of Angular Manifold", x = "w", 
                      size = 1.1)        
    if (target == "cdf")
      gg <- gg + labs(y = "Cross Section of Angular Manifold", x = "w", 
                      size = 1.1)  
    plot(gg)
    
  } else if (target == "Pickands") {
    frame <- data.frame(grid = c(0, fit$grid, 1), trajhat = c(1, fit$h, 1))
    ggplot(data = frame, aes(x = grid, y = trajhat)) +
      geom_line(aes(x = grid, y = trajhat), size = 1.1, col = "blue") +
      theme_classic() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 10)) +            
      coord_fixed(xlim = c(0, 1), ylim = c(.5, 1)) +
      labs(y = "Pickands Dependence Function", x = "w", size = 1.1)
  }  
}

