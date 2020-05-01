##  ========================================================================  ##
##  Miguel de Carvalho                                                        ##
##  Copyright (C) 2018                                                        ##
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

cmodes <- function(Y, thresholds = apply(Y[,-1], 2, quantile, probs = 0.95), 
                   nu = 100, ...)
    UseMethod("cmodes")

cmodes.default <- function(Y, thresholds = apply(Y[,-1], 2, quantile,
                           probs = 0.95), nu = 100, ...) {

    y <- Y[, -1]
    T <- dim(y)[1]
    K <- dim(y)[2]

    w <- list()
    k <- as.numeric()
    label <- as.numeric()
    v <- as.numeric()
    d <- as.numeric()
    for(j in 1:K) {
        w[[j]] <- which(y[, j] > thresholds[j]) / T
        k[j] <- length(w[[j]])
    }
    
    n <- sum(k)
    c <- list()
    istar <- array(NA, K)
    
    for(j in 1:K) {
        c[[j]] <- density(w[[j]], n = T, ...)
        istar[j] <- which(c[[j]]$y == max(c[[j]]$y))
    }
    
    m <- function(w)
      mean(dbeta(w, nu * istar / T, nu * (1 - istar / T)))
    M <- Vectorize(m)
    ## Organize and return outputs    
    outputs <- list(c = c, M = M, istar = istar, Y = Y, k = k)
    class(outputs) <- "cmodes"
    return(outputs)
}

plot.cmodes <- function(x, ylim = NULL, original = TRUE, sced = FALSE, ...) {
    T <- dim(x$Y)[1]
    K <- dim(x$Y)[2] - 1
    if (original == TRUE) {
        plot(x$Y[, 1], x$M(1:T / T), xlab = "Time (in years)", ylab = "Mode Mass Function", 
             type = "l", lwd = 2, ...)
        points(x$Y[x$istar, 1], rep(0, K), pch = 18)
        if(sced == TRUE) {
            for(j in 1:K) 
              lines(x$Y[, 1], x$c[[j]]$y, col = "grey")
            lines(x$Y[, 1], x$M(1:T / T), lwd = 2)
            }
    }
    if (original == FALSE) {
        par(pty = "s")
        plot(1:T / T, x$M(1:T / T), xlab = "Time", ylab = "Mode Mass Function", 
             type = "l", lwd = 2, ...)
        points(x$istar / T, rep(0, K), pch = 18)
        if(sced == TRUE) {
        for(j in 1:K)
          lines(x$c[[j]]$x, x$c[[j]]$y, col = "grey")
        lines(1:T / T, x$M(1:T / T), lwd = 2)
        }
    }
}