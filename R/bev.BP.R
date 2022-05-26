##  ========================================================================  ##
##  Miguel de Carvalho                                                        ##
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

bev.BP <- function(Y, tau = 0.95, c = 0.0001, m = 'maximize', T = 1000,
                   burn = 100, grid = seq(0.01, 0.99, length = 2^8), raw = TRUE,
                   target = "pdf")
    UseMethod("bev.BP")

bev.BP.default <- function(Y, tau = 0.95, c = 0.0001, m = 'maximize', T = 1000,
                           burn = 100, grid = seq(0.01, 0.99, length = 2^8), 
                           raw = TRUE, target = "pdf") {
    ## Basic input validation
    if (tau >= 1)
        stop('tau must be smaller than 1')
    if (2 * burn > T)
        stop('T must be larger than 2 * burn for adaptive MH')
    if (dim(Y)[2] > 2)
        stop('at the moment this function only applies to bivariate extremes')
    if (target != "pdf" & target != "cdf" & target != "Pickands" &
        target != "joint cdf")
        stop('target must be equal to "pdf", "cdf", "Pickands", or "joint cdf"')
    
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
    w <- cbind(w, 1 - w)
    k <- dim(w)[1]
    lgrid <- length(grid)

    p <- 2
    ## if (m == 'changepoint') {
    ##     m <- mable::mable(w[, 1], M = c(2, k), interval = c(0, 1))$m
    ##     J <- p + 2
    ##     while (m != choose(J - 1, p - 1)) 
    ##         J <- J + 1            
    ## }
    if (m == 'maximize') {
        ## Maximize J s.t. #basis functions <= #observations
        J <- p + 2
        m <- choose(J - 1, p - 1)
        while (m <= k) {
            J <- J + 1;
            m <- choose(J - 1, p - 1);
        }
        J <- J - 1
        m <- choose(J - 1, p - 1)
    } else if (is.numeric(m)) {
        J <- p + 2
        while (m != choose(J - 1, p - 1)) 
            J <- J + 1                                
    }
   
    ## Define R wrapper function to call the C function BERN
    bern_<- function(w, k, c, p, m, J, T, burn, grid, lgrid, target) {
        if (target == "pdf" || target == "regression")
            tar = 1
        if (target == "cdf")
            tar = 2
        
        .Call("BERN", as.double(w), as.integer(k), as.double(c), as.integer(p), 
              as.integer(m), as.integer(J), as.integer(T), as.integer(burn), 
              as.double(grid), as.integer(lgrid), as.integer(tar), 
              PACKAGE = "extremis")
    }
    
    Call <- bern_(w, k, c, p, m, J, T, burn, grid, lgrid, target)
    
    # 1 <= T <= burn: use static M-H update
    # burn < T <= 2 * burn: start building adaptive variance
    # 2 * burn < T: use adaptive variance and harvest iterates    
    PI <- Call$weights[(2 * burn + 1):T, ]
    PIhat <- colMeans(PI) 
    traj <- Call$traj[(2 * burn + 1):T, ]
    trajhat <- apply(traj, 2, mean) 
    
    ## Organize and return outputs    
    outputs <- list(PI = PI, PIhat = PIhat, indices = Call$indices, J = J, 
                    m = m, traj = traj, trajhat = trajhat, grid = grid, w = w, 
                    target = target)
    class(outputs) <- "bev.BP"
    return(outputs)
}

plot.bev.BP <- function(fit, bands = TRUE, target = 1) {
    grid <- fit$grid
    traj <- fit$traj
    trajhat <- fit$trajhat
    target <- fit$target

    frame <- data.frame(grid = grid, trajhat = trajhat)
    gg <- ggplot(data = frame, aes(x = grid, y = trajhat)) +
        geom_line(aes(x = grid, y = trajhat), size = 1.1, col = "blue") +
        theme_classic() +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 10))
    gg <- gg  + geom_rug(data = data.frame(w = fit$w[, 1]), aes(x = w, y = 0),
                         alpha = 0.5, sides = 'b')

    if (target == "pdf")
        gg <- gg + labs(y = "Angular Density", x = "w", size = 1.1)
    if (target == "cdf")
        gg <- gg + labs(y = "Angular Distribution", x = "w", size = 1.1)
    
    if (bands == TRUE) {
        l.band <- apply(traj, 2, quantile, probs = 0.025)
        u.band <- apply(traj, 2, quantile, probs = 0.975)
        gg <- gg + geom_ribbon(data = data.frame(l.band, u.band),
                               aes(ymin = l.band, ymax = u.band),
                               linetype = 2, alpha = 0.1)
    }
    plot(gg)
}
