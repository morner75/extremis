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

plotFrechet <- function(Y, tau = 0.95, raw = TRUE, ...)
    UseMethod("plotFrechet")

plotFrechet.default <- function(Y, tau = 0.95, raw = TRUE, ...) {
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
    par(pty = 's')
    plot(x, y, pch = 16, log = "xy", xlim = c(min(x, y), max(x, y)),
         ylim = c(min(x, y), max(x, y)), ...)
    u <- quantile(x + y, tau)
    X <- seq(-10000, 10000, .1)
    Y <- u - X
    points(X, Y, type = "l", lwd = 3, col = "grey")
    par(pty = 's')
}
