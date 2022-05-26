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

bev.kernel <- function(Y, tau = 0.95, nu, method = "euclidean",
                       grid = seq(0.01, 0.99, length = 2^8), raw = TRUE,
                       target = "pdf")
    UseMethod("bev.kernel")

bev.kernel.default <- function(Y, tau = 0.95, nu, method = "euclidean",
                               grid = seq(0.01, 0.99, length = 2^8),
                               raw = TRUE, target = "pdf") {
    ## Basic input validation
    if (tau >= 1)
        stop('tau must be smaller than 1')
    if(nu < 0)
        stop('nu must be positive')       
    if (dim(Y)[2] > 2)
        stop('this function only applies to bivariate extremes')
    if (target != "pdf" & target != "cdf" & target != "Pickands")
        stop('target must be equal to "pdf", "cdf", or "Pickands"')
    
    ## Convert data to common margins, if raw == TRUE
    if (raw == TRUE) {
        n <- dim(Y)[1]
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
    ## Compute angular distribution function
    k <- length(w)
    if (method == "euclidean") {
        v <- var(w) * ((k - 1) / k)
        p <- 1 / k * (1 - (mean(w) - 1/2) * v^{-1} * (w - mean(w)))
    }
    if (method == "empirical") {
        lambda <- el.test(w, 1 / 2)$lambda
        p <- 1 / k * (1 + (w - 1/2) * lambda)^-1
    }
    if (target == "pdf") {
        smooth <- function(W)
            sum(p * dbeta(W, w * nu, (1 - w) * nu))        
        estimate <- sapply(grid, smooth)
    }
    if (target == "cdf") {
        smooth <- function(W)
            sum(p * pbeta(W, w * nu, (1 - w) * nu))
        estimate <- sapply(grid, smooth)
    }
    if (target == "Pickands") {
        smooth <- function(W) {
            ptilde <- numeric(k)
            integrand <- function(W)
                pbeta(W, nu * w[i], (1 - w[i]) * nu)
            for(i in 1:k) 
                ptilde[i] <- p[i] * as.numeric(integrate(integrand, 0, W)[1])            
            1 - W + 2 * sum(ptilde)
        }
        smooth <- Vectorize(smooth)
        estimate <- sapply(grid, smooth)
    }
    if (target == "joint cdf") {
        ibeta <- function(x, a, b)
            pbeta(x, a, b) * beta(a, b) 
        smooth <- function(x, y) {
            ptilde <- numeric(k)
            ## integrand <- function(W)
            ##     max(W / x, (1 - W) / y) * dbeta(W, nu * w[i], (1 - w[i]) * nu)
            for(i in 1:k)
                ptilde[i] <- p[i] / x *
                    (w[i] - ibeta(x / (x + y), nu * w[i] + 1, (1 - w[i]) * nu) /
                     beta(nu * w[i], (1 - w[i]) * nu)) + 
                    p[i] / x * (ibeta(x / (x + y), nu * w[i], (1 - w[i]) * nu) -
                                ibeta(x / (x + y), nu * w[i] + 1, (1 - w[i]) * nu) / 
                                beta(nu * w[i], (1 - w[i]) * nu))
            ## ptilde[i] <- p[i] * as.numeric(integrate(integrand, 0, 1)[1])
            exp(-2 / k * sum(ptilde))      
        }
        max.x <- max.y <- 10
        xx <- seq(1, max.x, 1)
        yy <- seq(1, max.y, 1)
        estimate <- matrix(0, max.x, max.y)
        for (i in 1:max.x)        
            for (j in 1:max.y)       
                estimate[i, j] <- smooth(xx[i], yy[j])
    }    
    
    ## Organize and return outputs
    outputs <- list(estimate = estimate, p = p, grid = grid, w = w, nu = nu, 
                     target = target)
    ## outputs <- list(estimate = estimate, p = p, w = w, nu = nu, grid = grid,
    ##                 Y = Y, Frechet = cbind(x, y), target = target)
    class(outputs) <- "bev.kernel"
    return(outputs)
}

plot.bev.kernel <- function(fit, xlab = "w", ...) {
    grid <- fit$grid
    trajhat <- fit$estimate
    target <- fit$target

    if (target == "pdf" | target == "cdf") {
        frame <- data.frame(grid = grid, trajhat = trajhat)
        gg <- ggplot(data = frame, aes(x = grid, y = trajhat)) +
            geom_line(aes(x = grid, y = trajhat), size = 1.1, col = "blue") +
            theme_classic() +
            theme(axis.text = element_text(size = 10),
                  axis.title = element_text(size = 10))
        gg <- gg  + geom_rug(data = data.frame(w = fit$w), aes(x = w, y = 0),
                             alpha = 0.5, sides = 'b')    
        if (target == "pdf")
            gg <- gg + labs(y = "Angular Density", x = "w", size = 1.1)
        if (target == "cdf")
            gg <- gg + labs(y = "Angular Distribution", x = "w", size = 1.1)        
        plot(gg)
    } else if (target == "Pickands") {
        frame <- data.frame(grid = c(0, fit$grid, 1), trajhat = c(1, fit$estimate, 1))
        ggplot(data = frame, aes(x = grid, y = trajhat)) +
            geom_line(aes(x = grid, y = trajhat), size = 1.1, col = "blue") +
            theme_classic() +
            theme(axis.text = element_text(size = 10),
                  axis.title = element_text(size = 10)) +            
            coord_fixed(xlim = c(0, 1), ylim = c(0.5, 1)) +
            labs(y = "Pickands Dependence Function", x = "w", size = 1.1)         
        ## lines(c(0,1), c(1, 1), lty = 3)
        ## lines(c(0,.5), c(1,.5), lty = 3)
        ## lines(c(0.5, 1), c(.5, 1), lty = 3)
    }
}
