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

angcdf <- function(Y, tau = 0.95, method = "euclidean", raw = TRUE)
    UseMethod("angcdf")

angcdf.default <- function(Y, tau = 0.95, method = "euclidean", raw = TRUE) {
    ## Basic input validation
    if (tau >= 1)
        stop('tau must be smaller than 1')
    if (dim(Y)[2] > 2)
        stop('this function only applies to bivariate extremes')
    
    ## Convert data to common margins, if raw == TRUE
    if(raw == TRUE) {
        n <- dim(Y)[1]
        FX <- ecdf(Y[, 1])
        FY <- ecdf(Y[, 2])
        x <- -1 / log(n / (n + 1) * FX(Y[, 1]))
        y <- -1 / log(n / (n + 1) * FY(Y[, 2]))
    }
    if(raw == FALSE) {
        x <- Y[, 1]
        y <- Y[, 2]
    }
    ## Compute angles and threshold radius
    w <- x / (x + y)
    u <- quantile(x + y, tau)    
    w <- w[which(x + y > u)]
    ## Compute angular distribution function
    w <- sort(w) 
    k <- length(w)
    if(method == "euclidean") {
        v <- var(w) * (k / (k - 1))
        p <- 1 / k * (1 - (mean(w) - 1 / 2) * v^{-1} * (w - 1 / 2))
    }
    if(method == "empirical") {
        delta <- el.test(w, 1 / 2)$lambda
        p <- 1 / k * (1 + (w - 1 / 2) * delta)^-1
    }
    H <- approxfun(w, cumsum(p), method = "constant", yleft = 0,
                   yright = 1, f = 0, ties = "ordered")    
    ## Organize and return outputs    
    outputs <- list(H = H, w = w, Y = Y)
    class(outputs) <- "angcdf"
    return(outputs)
}

plot.angcdf <- function(x, ylim = NULL, original = TRUE, main = "", ...) 
    plot(x$H, xlab = "w", ylab = "Angular Distribution Function",
         main = "", type = "S", ...)

