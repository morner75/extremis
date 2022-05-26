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

scedasis.kernel <- function(XY, tau = 0.95, raw = TRUE, structure = "min", ...)
    UseMethod("scedasis.kernel")

scedasis.kernel.default <- function(XY, tau = 0.95, raw = TRUE,
                                    structure = "min", ...) {

    if ((is.matrix(XY) | is.data.frame(XY)) == FALSE)
        stop('data needs to be a matrix or data frame')
    dim <- ncol(XY)
    if (dim == 2) {
        y <- XY[,2]
    } else {
        ## Convert data to common margins, if raw == TRUE
        if(raw == TRUE) {
            n <- dim(XY)[1]
            aux <- function(c)
                function(c) -1 / log(length(c) / (length(c) + 1) * ecdf(c)(c))
            FZ <- apply(XY[,-1], 2, aux)
        }
        if(raw == FALSE) 
            FZ<-XY[,-1]
        if(structure == "min")
            y <- apply(FZ,1,min)
        if(structure == "max")
            y <- apply(FZ,1,max)
        if(structure == "sum")
            y <- apply(FZ,1,sum)
    }
    threshold <- quantile(y, tau)
    ## Basic input validation
    if (as.numeric(threshold) >= max(y))
        stop('threshold cannot exceed max(y)')
    
    ## Initialize variables
    T <- length(y)
    w <- which(y > threshold) / T
    k <- length(w) 
    c <- density(w, n = T, from = 1 / T, to = 1, ...)
    
    ## Organize and return outputs    
    outputs <- list(c = c, w = w, k = k, T = T, XY = XY)
    class(outputs) <- "scedasis.kernel"
    return(outputs)
}

plot.scedasis.kernel <- function(x, rugrep = TRUE, original = TRUE,
                                 main = "", ...) {
    dim <- ncol(x$XY)
    if(dim == 2){
        labs<-"Scedasis Density"
    } else {
        labs<-"Structure Scedasis Density"
    }
    if(original == TRUE) {
        plot(x$XY[, 1], x$c$y, xlab = "Time", ylab = labs, main = "",
             type = "S", ...)
        if(rugrep == TRUE)
            rug(x$XY[x$w * x$T, 1])
    }
    if(original == FALSE) {
        par(pty = "s")
        plot(x$c, xlab = "w", ylab = labs, main = "", ...)
        if (rugrep == TRUE)
            rug(x$w)
    }
}
