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

cdensity <- function(Y, threshold = quantile(Y[, 2], 0.95), ...)
    UseMethod("cdensity")

cdensity.default <- function(Y, threshold = quantile(Y[, 2], 0.95), ...) {
    y <- Y[, 2]
    ## Basic input validation
    if (as.numeric(threshold) >= max(y))
        stop('threshold cannot exceed max(y)')

    ## Initialize variables
    T <- length(y)
    w <- which(y > threshold) / T
    k <- length(w) 
    c <- density(w, n = T, from = 1 / T, to = 1, ...)
    
    ## Organize and return outputs    
    outputs <- list(c = c, w = w, k = k, T = T, Y = Y)
    class(outputs) <- "cdensity"
    return(outputs)
}

plot.cdensity <- function(x, rugrep = TRUE,
                          original = TRUE, main = "", ...) {
    if(original == TRUE) {
        plot(x$Y[, 1], x$c$y, xlab = "Time", ylab = "Scedasis Density", 
             main = "", type = "S", ...)
        if(rugrep == TRUE)
            rug(x$Y[x$w * x$T, 1])
    }
    if(original == FALSE) {
        par(pty = "s")
        plot(x$c, xlab = "w", ylab = "Scedasis Density", 
             main = "", ...)
        if (rugrep == TRUE)
          rug(x$w)
    }
}
