##  ========================================================================  ##
##  Miguel de Carvalho and Junho Lee                                          ##
##  Copyright (C) 2022                                                        ##
##  ------------------------------------------------------------------------  ##
##  This program is free software; you can redistribute it and/or modify      ##
##  it under the terms of the GNU General Public License as published by      ##
##  the Free Software Foundation; either version 2 of the License, or         ##
##  (at your option) any later version.                                       ##
##                                                                            ##
##  This program is distributed in the hope that it will be useful,           ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of            ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             ##
##  GNU General Public License for more details.                              ##
##                                                                            ##
##  You should have received a copy of the GNU General Public License         ##
##  along with this program; if not, a copy is available at                   ##
##  http://www.r-project.org/Licenses/                                        ##
##  ========================================================================  ##

chit <- function(XY, tau = 0.95, prior = list(a = 0.001, b = 0.001),
                 T = 10000, burn = 5000, link = 'probit', knot = 20, dif = 2,
                 raw = TRUE, silent = FALSE)
  UseMethod("chit")



chit.default <- function(XY, tau = 0.95, prior = list(a = 0.001, b = 0.001),
                         T = 10000, burn = 5000, link = 'probit', knot = 20,
                         dif = 2, raw = TRUE, silent = FALSE) {
  
  ## Basic input validation
  if (tau >= 1)
    stop('tau must be smaller than 1')
  if (dim(XY)[2] > 3)
    stop('this function only applies to bivariate extremes')
  
  ## Convert data to common margins, if raw == TRUE
  time <- XY[, 1]    
  if (raw == TRUE) {
    n <- dim(XY)[1]
    FX <- ecdf(XY[, 2])
    FY <- ecdf(XY[, 3])
    x <- -1 / log(n / (n + 1) * FX(XY[, 2]))
    y <- -1 / log(n / (n + 1) * FY(XY[, 3]))
  }
  if (raw == FALSE) {
    x <- XY[, 2]
    y <- XY[, 3]
  }
  u <- quantile(y, tau)
  time <- time[which(y > u)]
  IE <- ifelse(x[which(y > u)] > u, 1, 0)
  k <- length(IE)        
  
  ## PRELIMINARIES ON P-SPLINES
  ## penalty matrix function
  Kdmat <- function(j, d) {
    Dd <- cbind(0, diag(j - 1)) - cbind(diag(j - 1), 0)
    if (d == 1) {
      return(t(Dd) %*% Dd)
    } else { 
      for(i in 2:d) {
        Dd <- (cbind(diag(j - d), 0) - cbind(0, diag(j - d))) %*% Dd}
    }
    return(t(Dd) %*% Dd)
  }    
  ## B-spline basis matrix and evaluate penalty matrix    
  Bs <- bs(time, knots = seq(min(time), max(time), length = knot),
           degree = 3, intercept = TRUE)   
  K <- Kdmat(ncol(Bs), dif)
  output <- matrix(nrow = T, ncol = ncol(Bs) + 1)
  nbs <- ncol(Bs)
  
  ## PROBIT LINK: SET INITIAL VALUES AND CALL GIBBS SAMPLER    
  if (link == "probit") {
    t2 <- 1
    beta <- rnorm(nbs)    
    Sigmaf <- function(t2) solve((t(Bs) %*% Bs) + K / t2, tol = 2.220446e-100)        
    if(!silent){
      cat("MCMC iterations:\n")
      cat("================\n") 
    }       
    for (i in 1:T) {
      ## update latent indicators (li)
      a <- ifelse(IE == 0, -Inf, 0)
      b <- ifelse(IE == 1, Inf, 0)
      li <- truncnorm::rtruncnorm(length(IE), a, b, Bs %*% beta, sd = 1)        
      ## update beta
      beta <- as.vector(rmvnorm(1, mean = Sigmaf(t2) %*% (t(Bs) %*% li),sigma = Sigmaf(t2)))

      ## update t2
      t2 <- 1 / rgamma(1, prior$a + (length(beta) - dif) / 2,
                       prior$b + t(beta) %*% (K %*% beta) / 2)
      
      output[i, 1:nbs] <- beta
      output[i, nbs+1] <- t2
      if (!silent & i%%(T %/% 5) == 0)
        cat(i,"\n")
    }
  }
  ## PROBIT AND CLOGOG LINKS: SET INITIAL VALUES AND CALL METROPOLIS-HASTINGS
  if (link == "logit" | link == "cloglog") {
    
    linkf <- make.link(link=link)
    
    ## auxiliary functions for weight matrix and working response
    lik <- function(Bsb)
      sum(dbinom(IE, 1, prob = linkf$linkinv(Bsb),log = TRUE))
    W <- function(Bsb) 
      diag(as.vector(linkf$mu.eta(Bsb)/linkf$linkinv(Bsb)*(1-linkf$linkinv(Bsb))))
    ytilde <- function(Bsb)  Bsb+ (IE - linkf$linkinv(Bsb))/linkf$mu.eta(Bsb)
    
    ## auxiliary functions for mean and covariance matrix
    meanc <- function(co, w, Bsb)
      as.vector(co %*% (t(Bs) %*% (w %*% ytilde(Bsb))))
    covc <- function(w) 
      solve(t(Bs) %*% w %*% Bs + K / t2, tol = 2.220446e-100)
      
    
    ## SET INITIAL VALUES
    ## initialize at the posterior mode using conjugate gradient
    beta <- rep(0, nbs)
    ran <- (length(beta) - dif)
    init <- function(x) {
      beta <- x[1:nbs]
      Bsb <- Bs%*%beta
      sum(dbinom(IE, 1,
                 prob = linkf$linkinv(Bsb),
                 log = TRUE))  +
        as.vector(-ran / 2  -
                    t(beta) %*% K %*% beta / 2)
    }
    init_derv <- function(x) {
      beta <- x[1:nbs]
      Bsb <- Bs%*%beta
      t(Bs)%*%W(Bsb)%*%
                 (IE-linkf$linkinv(Bsb))- K%*%beta/2
    }
    cat("Initializing...\n")
    x <- optim(beta, init,init_derv,
               control = list(fnscale = -1000, reltol = 10^-20,
                              maxit = 10^9), method = "BFGS")
    
    if (x$convergence == 0)
      cat("(converged)\n")
    beta <- x$par[1:nbs]
    t2 <- 1#exp(x$par[nbs + 1])                


    ## START METROPOLIS-HASTINGS
    cat("MCMC iterations:\n")
    cat("================\n")
    for (i in 1:T) {
      ## update beta
      Bsb <- Bs %*% beta
      w <- W(Bsb)
      if (is.nan(det(t(Bs) %*% w %*% Bs + K / t2)) == FALSE &
          det(t(Bs) %*% w %*% Bs + K / t2) != 0) {            
        cbetap <- covc(w)
        ## if not numerically symmetric, convert it to symmetric
        if (base::isSymmetric(cbetap) == FALSE)
          cbetap[lower.tri(cbetap)] <- t(cbetap)[lower.tri(cbetap)] 
        
        cbeta <- covc(w)
        mbeta <- meanc(cbeta, w, Bsb)
        
        betap <- as.vector(mvtnorm::rmvnorm(1, mean = mbeta, sigma = cbeta))
        Bsbp <- Bs %*% betap
        w <- W(Bsbp)            
        cbetap <- covc(w)
        mbetap <- meanc(cbetap, w, Bsbp)
        
        like1c <- lik(Bsb) 
        like1p <- lik(Bsbp)
        like2c <- -ran / 2 * log(t2) - 
          t(beta) %*% K %*% beta / (2 * t2)
        like2p <- -ran / 2 * log(t2) - 
          t(betap) %*% K %*% betap / (2 * t2)
        evalp <- mvtnorm::dmvnorm(betap, mean = mbeta, sigma = cbeta,
                                  log = TRUE)
        evalc <- mvtnorm::dmvnorm(beta, mean = mbetap, sigma = cbetap,
                                  log = TRUE)
        alpha <- like1p + like2p - evalp  - (like1c + like2c - evalc) 
        
        if (is.nan(alpha) == FALSE) 
          if (alpha > log(runif(1)))
            beta <- betap
      }
      ## update t2
      t2 <- 1 / rgamma(1, prior$a + (length(beta) - dif) / 2,
                       prior$b + t(beta) %*%
                         Kdmat(length(beta), dif) %*% beta / 2)
      output[i, ] <- c(beta, t2)
      if (i%%(T %/% 5) == 0)
        cat(i,"\n")            
    }    
  }
  
  beta.post <-  output[(burn + 1):T, 1:ncol(Bs)]
  t2.post <- output[(burn + 1):T, ncol(Bs) + 1]
  if (link == 'probit')
    traj <- pnorm(beta.post %*% t(Bs))
  if (link == 'logit')
    traj <- 1 / (1 + exp(-beta.post %*% t(Bs)))
  if (link == 'cloglog')
    traj <- 1 - exp(-exp(beta.post %*% t(Bs)))
  if(!silent) cat("DONE\n")
  trajhat <- colMeans(traj)
  outputs <- list(IE = IE, tex = time, traj = traj, trajhat = trajhat, 
                  beta = beta.post, t2 = t2.post,
                  afun = approxfun(time, trajhat))
  class(outputs) <- "chit"
  return(outputs)
}


plot.chit <- function(fit, bands = TRUE) {
    frame <- data.frame(x = fit$tex, y = fit$IE, trajhat = fit$trajhat)    
    p <- ggplot(data = frame, aes(x = x, y = y)) +
        geom_point(aes(x = x, y = y), shape = 16, alpha = 0.5) +
        geom_line(aes(x = x, trajhat), size = 1.1, col = "blue") +
        xlab("Time") +
        ylab(expression(chi[t])) +
        ## theme_classic() +
        theme_bw() + 
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 10))    
    if (bands == TRUE) {
        l.band <- apply(fit$traj, 2, quantile, probs = 0.025)
        u.band <- apply(fit$traj, 2, quantile, probs = 0.975)        
        p <- p +
            geom_ribbon(aes(ymin = l.band, ymax = u.band),
                        linetype = 2, alpha = 0.1)
    }
    plot(p)    
}


