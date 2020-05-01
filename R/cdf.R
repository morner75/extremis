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

cdf <- function(Y, threshold = quantile(Y[, 2], 0.95))
    UseMethod("cdf")

cdf.default <- function(Y, threshold = quantile(Y[, 2], 0.95)) {
    y <- Y[, 2]
    ## Basic input validation
    if (as.numeric(threshold) >= max(y))
        stop('threshold cannot exceed max(y)')

    ## Initialize variables
    T <- length(y)
    I <- 1:T
    w <- which(y > threshold) / T
    k <- length(w) 
    C <- ecdf(w)
    
    ## Organize and return outputs    
    outputs <- list(C = C, w = w, Y = Y)
    class(outputs) <- "cdf"
    return(outputs)
}

plot.cdf <- function(x, uniform = TRUE, ylim = NULL,
                     original = TRUE, main = "", ...) {
    if(original == TRUE) {
        T <- dim(x$Y)[1]
        I <- 1:T
        plot(x$Y[, 1], x$C(I / T), xlab = "Time",
             ylab = "Scedasis Distribution Function", 
             main = "", type = "S", ...)
        if(uniform == TRUE)
            lines(x$Y[, 1], seq(0, 1, length = T), lty = 2)
    }
    if(original == FALSE) {
        par(pty = "s")
        plot(x$C, xlab = "w", ylab = "Scedasis Distribution Function", 
             main = "", ...)
        if(uniform == TRUE)
            abline(0, 1, lty = 2)        
    }
}
