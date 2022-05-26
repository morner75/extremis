##  ========================================================================  ##
##  Miguel de Carvalho and Alina Kumukova                                     ##
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

regmani.BP <- function(Y, tau = 0.95, c = 0.0001, m = 'maximize', T = 1000,
                       burn = 100, raw = TRUE)
    UseMethod("regmani.BP")

regmani.BP.default <- function(Y, tau = 0.95, c = 0.0001, m = 'maximize',
                               T = 1000, burn = 100, raw = TRUE) {
    ## Input validation
    if (tau >= 1)
        stop('tau must be smaller than 1')
    if (2 * burn > T)
        stop('T must be larger than 2 * burn for adaptive MH')
    if (dim(Y)[2] > 2)
        stop('at the moment this function only applies to bivariate extremes')
    
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
    ## Store data in unit Frechet margins
    Y.unitFrechet <- cbind(x,y)
    ## Compute angles and threshold radius
    w <- x / (x + y)
    u <- quantile(x + y, tau)
    exceedances <- cbind(x[x + y > u], y[x + y > u])
    w <- w[which(x + y > u)]
    w <- cbind(w, 1 - w)
    k <- dim(w)[1]

    p <- 2
    if (m == 'changepoint') {
        m <- mable::mable(w[, 1], M = c(2, k), interval = c(0, 1))$m
        J <- p + 2
        while (m != choose(J - 1, p - 1)) 
            J <- J + 1            
    } else if (m == 'maximize') {
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
        .Call("BERN", as.double(w), as.integer(k), as.double(c), as.integer(p), 
              as.integer(m), as.integer(J), as.integer(T), as.integer(burn), 
              as.double(grid), as.integer(lgrid), as.integer(target), 
              PACKAGE = "extremis")
    }
    target <- 0
    grid <- seq(0.01, 0.99, length = 2^8)
    lgrid <- length(grid)
    Call <- bern_(w, k, c, p, m, J, T, burn, grid, lgrid, target)
    
    # 1 <= T <= burn: use static M-H update
    # burn < T <= 2 * burn: start building adaptive variance
    # 2 * burn < T: use adaptive variance and harvest iterates    
    PI <- Call$weights[(2 * burn + 1):T, ]
    allPI <- Call$weights
    PIhat <- apply(PI, 2, mean)
    
    ## Organize and return outputs    
    outputs <- list(PI = PI, PIhat = PIhat, indices = Call$indices, J = J, 
                    m = m, w = w, Y  = Y, exceedances = exceedances, Y.unitFrechet = Y.unitFrechet)
    class(outputs) <- "regmani.BP"
    return(outputs)
}

plot.regmani.BP <- function(fit, tau, q = NULL, pred = NULL, bands = FALSE, points = TRUE, 
                            original = FALSE, qqplot = FALSE, cex = 1,
                            which = 'manifold', ...) {
    ## Regression manifold workhorse functions
    Gcond2D <- function(x, y, w, q) {
        ## conditional distribution estimate for Y | X = x  
        a1 <- w * fit$indices[, 1]
        a2 <- w * fit$indices[, 2]
        b1 <- 1 - pbeta(x / (x + y), fit$indices[, 1] + 1, fit$indices[, 2])
        b2 <- pbeta(x / (x + y), fit$indices[, 1], fit$indices[, 2] + 1)
        return(2 / fit$J * exp(-2 / fit$J * (x^{-1} * (a1%*%b1) + y^{-1} *
                                             (a2%*%b2)) + x^{-1}) * (a1%*%b1) - q)
    }
    
    reg <- function(xq, PI) 
        uniroot(Gcond2D, c(0.0005, 10^12), x = xq[1],
                w = PI, q = xq[2], maxiter = 10000, check.conv = TRUE, tol = 1e-15)$root
    
    if (which == 'manifold') {
        ## Plot the fitted regression manifold
        if (is.null(q) == TRUE & is.null(pred) == TRUE) {
            l <- 30
            if (original == FALSE) {
                ## grid over which the regression manifold is evaluated           
                gr <- expand.grid(x = seq(0.005, 20, length = l),
                                  q = seq(0.01, 0.99, length = l))                
                RM <- data.frame(x = gr[, 1], y = gr[, 2], 
                                 z = apply(gr, 1, reg, PI = fit$PIhat))                
                mg <- data.frame('x2' = fit$Y.unitFrechet[, 1], 'z2' = fit$Y.unitFrechet[, 2])
                mg$y2 <- rep(0, nrow(mg))                
                xlim <- c(0, 20)
            } else if (original == TRUE) {
                gr <- expand.grid(x = seq(0.01, max(fit$Y.unitFrechet[, 1]), length = l),
                                  q = seq(0.01, 0.99, length = l))                
                RM <- data.frame(x = gr[, 1], y = gr[, 2], 
                                 z = apply(gr, 1, reg, PI = fit$PIhat))                
                mg <- data.frame('x2' = fit$Y[, 1], 'z2' = fit$Y[, 2])
                mg$y2 <- rep(0, nrow(mg))
                RM$x <- quantile(ecdf(fit$Y[, 1]), prob = exp(-1 / RM$x))
                RM$z <- quantile(ecdf(fit$Y[, 2]), prob = exp(-1 / RM$z))
                xlim <- range(fit$Y[, 1])
                zlim <- range(fit$Y[, 2])
            }

            blues_pal <- colorRampPalette(c("darkblue", 'blue', 'cyan',
                                            'green4', 'yellow'), bias = 1.5)(500)
            blues_pal_transp <- transparent(blues_pal, trans.val = 0.4)
            
            mypanel <- function(x, y, z, x2, y2, z2,...) {
                panel.wireframe(x, y, z,...)
                panel.cloud(x2, y2, z2, pch = 16, col = 'darkblue', ...)
            }
            
            wireframe(z ~ x * y, data = RM,
                      scales = list(arrows = FALSE, col = "black", cex = cex),
                      xlab = list("x", cex = cex, rot = -30),
                      ylab = list("q", cex = cex),
                      zlab = list(latex2exp::TeX('$y_{q | x}$'), rot = 90, cex = cex),
                      xlim = xlim,                      
                      ylim = c(0, 1),
                      par.settings = list(axis.line = list(col = "transparent")),
                      drape = T, colorkey = F, par.box = list(col = 1, lty = 1),
                      col.regions = blues_pal_transp,
                      screen = list(z = 120, x = -80, y = 5),
                      panel.3d.wireframe =
                          function(x, y, z,
                                   xlim, ylim, zlim,
                                   xlim.scaled, ylim.scaled, zlim.scaled,
                                   pts,
                                   ...) {
                              panel.3dwire(x = x, y = y, z = z,
                                           xlim = xlim,
                                           ylim = ylim,
                                           zlim = zlim,
                                           xlim.scaled = xlim.scaled,
                                           ylim.scaled = ylim.scaled,
                                           zlim.scaled = zlim.scaled,
                                           ...)
                              xx <-
                                  xlim.scaled[1] + diff(xlim.scaled) *
                                  (mg$x2 - xlim[1]) / diff(xlim)
                              yy <-
                                  ylim.scaled[1] + diff(ylim.scaled) *
                                  (mg$y2 - ylim[1]) / diff(ylim)
                              zz <-
                                  zlim.scaled[1] + diff(zlim.scaled) *
                                  (mg$z2 - zlim[1]) / diff(zlim)
                              panel.3dscatter(x = xx,
                                              y = yy,
                                              z = zz,
                                              xlim = xlim,
                                              ylim = ylim,
                                              zlim = zlim,
                                              xlim.scaled = xlim.scaled,
                                              ylim.scaled = ylim.scaled,
                                              zlim.scaled = zlim.scaled,
                                              pch = 16,col = 'darkblue',
                                              ...)
                          }, ...)        
        } else if (is.numeric(q) == TRUE) {
            ## Basic input validation
            if (sum(q >= 1) > 0  || sum(q < 0) > 0)
                stop('q must be inside the unit interval')
            
            iterations <- dim(fit$PI)[1]
            xmin <- 0.005
            
            if (original == FALSE) {
                
                xmax <- 20       
                grp <- expand.grid(x = seq(xmin, xmax, length = 50), q = q)
                csp <- data.frame(x = grp[, 1], q = grp[, 2],
                                  Zest = apply(grp, 1, reg, PI = fit$PIhat))
            
                if (bands == TRUE) {
                    cat('Learning...\n')            
                    cores <- detectCores()
                    cat('using ', cores, ' cores...\n')            
                    cl <- makeCluster(cores[1] - 1)
                    print(cl)
                    registerDoParallel(cl)            
                    cip <- foreach(i = 1:iterations, .combine = cbind) %dopar% {            
                        tempMatrix <- apply(grp, 1, reg, PI = fit$PI[i, ])                
                        tempMatrix
                    }
                    stopCluster(cl)
                    csp$lower <- apply(cip, 1, quantile, 0.025)
                    csp$upper <- apply(cip, 1, quantile, 0.975)
                    csp$q <- as.factor(csp$q)
                    cat('DONE \n')
                }
                
                mg <- data.frame('x2' = fit$Y.unitFrechet[, 1], 'z2' = fit$Y.unitFrechet[, 2])
                mg$y2 <- rep(0, nrow(mg))
                
                p <- ggplot(data = csp, aes(x = x, y = Zest, col = as.factor(q)))
                    
                if (points == TRUE)
                    p <- p +
                        geom_point(data = mg, aes(x = x2, y = z2), color = "gray57", pch = 16) 
                
                p <- p + 
                    geom_line(aes(y = Zest, group = as.factor(q)), size = 1.7) +
                    theme_classic() +
                    coord_fixed(xlim = c(0, 20), ylim = c(0, 20)) +
                    labs(y = "y", x = "x") +
                    guides(col = guide_legend(title = "q", order = 2), 
                           linetype = guide_legend(title="Regression line", order = 1)) + 
                    theme(axis.text = element_text(size = 15), 
                          axis.title = element_text(size = 15),
                          legend.text = element_text(size = 15),
                          legend.title = element_text(size = 15), 
                          legend.key.size = unit(2, "line"))
                
                if (bands == TRUE)
                    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper), 
                                         linetype = 0, alpha = 0.1)
                plot(p)
                
            } else if (original == TRUE) {
                
                xmax <- max(fit$Y.unitFrechet[,1])    
                grp <- expand.grid(x = seq(xmin, xmax, length = 75), q = q)
                csp <- data.frame(x = grp[, 1], q = grp[, 2],
                                  Zest = apply(grp, 1, reg, PI = fit$PIhat))
                
                if (bands == TRUE) {
                    cat('Learning...\n')            
                    cores <- detectCores()
                    cat('using ', cores, ' cores...\n')            
                    cl <- makeCluster(cores[1] - 1)
                    print(cl)
                    registerDoParallel(cl)            
                    cip <- foreach(i = 1:iterations, .combine = cbind) %dopar% {            
                        tempMatrix <- apply(grp, 1, reg, PI = fit$PI[i, ])                
                        tempMatrix
                    }
                    stopCluster(cl)
                    csp$Zest <- apply(cip, 1, mean)
                    for (i in 1:dim(cip)[2]) cip[,i]=quantile(ecdf(fit$Y[,2]),prob=exp(-1/cip[,i]))
                    csp$lower <- apply(cip, 1, quantile, 0.025)
                    csp$upper <- apply(cip, 1, quantile, 0.975)
                    csp$q <- as.factor(csp$q)
                    cat('DONE \n')
                }
                
                mg <- data.frame('x2' = fit$Y[, 1], 'z2' = fit$Y[, 2])
                mg$y2 <- rep(0, nrow(mg))
                
                p <- ggplot(data = csp, aes(x = quantile(ecdf(fit$Y[,1]),prob=exp(-1/x)), y = quantile(ecdf(fit$Y[,2]),prob=exp(-1/Zest)), col = as.factor(q)))
                
                if (points == TRUE)
                    p <- p +
                        geom_point(data = mg, aes(x = x2, y = z2), color = "gray57", pch = 16) +
                        geom_point(data = mg[rowSums(fit$Y.unitFrechet) > quantile(rowSums(fit$Y.unitFrechet),probs=tau),],
                                   aes(x = x2, y = z2), color = "black", pch = 16) 
                
                p <- p + 
                    geom_line(aes(y = quantile(ecdf(fit$Y[,2]),prob=exp(-1/Zest)), group = as.factor(q)), size = 1.7) +
                    theme_classic() +
                    coord_fixed() +
                    labs(y = "y", x = "x") +
                    guides(col = guide_legend(title = "q", order = 2), 
                           linetype = guide_legend(title="Regression line", order = 1)) + 
                    theme(axis.text = element_text(size = 15), 
                          axis.title = element_text(size = 15),
                          legend.text = element_text(size = 15),
                          legend.title = element_text(size = 15), 
                          legend.key.size = unit(2, "line"))
                
                if (bands == TRUE)
                    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper), 
                                         linetype = 0, alpha = 0.1)
                plot(p)
                
            }
            
        } else if (is.numeric(pred) == TRUE) {
            
            iterations <- dim(fit$PI)[1]
            qmin <- 0.01
            qmax <- 0.99        
            
            if (original == FALSE) {
                x <- pred
                grx <- expand.grid(x = x, q = seq(qmin, qmax, length = 50))
                csx <- data.frame(x = grx[, 1], q = grx[, 2],
                              Zest = apply(grx, 1, reg, PI = fit$PIhat))
            
                if (bands == TRUE) {
                    cat('Learning...\n')            
                    cores <- detectCores()
                    cat('using ', cores, ' cores...\n')            
                    cl <- makeCluster(cores[1] - 1)
                    print(cl)
                    registerDoParallel(cl)
                    cix <- foreach(i = 1:iterations, .combine=cbind) %dopar% {
                        tempMatrix <- apply(grx, 1, reg, PI = fit$PI[i, ])                           
                        tempMatrix           
                    }
                    stopCluster(cl)                
                    csx$lower <- apply(cix, 1, quantile, 0.025)
                    csx$upper <- apply(cix, 1, quantile, 0.975)
                    csx$x <- as.factor(csx$x)                
                    cat('DONE \n')
                }
                
                p <- ggplot(data = csx, aes(x = q, y = Zest, col = as.factor(x))) +
                    geom_line(aes(y = Zest, group = as.factor(x)), size = 1.7) +
                    theme_classic() +
                    labs(y = "Conditional quantile", x = "q") +
                    guides(linetype = guide_legend(title=" "), col = guide_legend(title = "x")) + 
                    theme(axis.text = element_text(size = 15),
                          axis.title = element_text(size = 15),
                          legend.text = element_text(size = 15),
                          legend.title = element_text(size = 15), legend.key.size = unit(2, "line"))
                plot(p)
                
                if (bands == TRUE)
                    p <- p + geom_ribbon(aes(ymin = lower, ymax = upper), 
                                         linetype = 0, alpha = 0.1)
                plot(p)
                
            } else if (original == TRUE){
                
                # convert to Frechet scale
                FX <- ecdf(fit$Y[, 1])
                x <- -1 / log(nrow(fit$Y) / (nrow(fit$Y) + 1) * FX(pred))
                grx <- expand.grid(x = x, q = seq(qmin, qmax, length = 50))
                csx <- data.frame(x = grx[, 1], q = grx[, 2],
                                  Zest = apply(grx, 1, reg, PI = fit$PIhat))
                
                csx$x <- expand.grid(x = pred, q = seq(qmin, qmax, length = 50))[,1]
                csx$Zest <- quantile(ecdf(fit$Y[,2]),prob=exp(-1/csx$Zest))
                
                p <- ggplot(data = csx, aes(x = q, y = Zest, col = as.factor(x))) +
                    theme_classic() +
                    coord_cartesian(xlim=c(0,1)) +
                    labs(y = "Conditional quantile", x = "q") +
                    guides(linetype = guide_legend(title=" "), col = guide_legend(title = "x")) + 
                    theme(axis.text = element_text(size = 15),
                          axis.title = element_text(size = 15),
                          legend.text = element_text(size = 15),
                          legend.title = element_text(size = 15), legend.key.size = unit(2, "line"))
                
                if (bands == TRUE) {
                    cat('Learning...\n')            
                    cores <- detectCores()
                    cat('using ', cores, ' cores...\n')            
                    cl <- makeCluster(cores[1] - 1)
                    print(cl)
                    registerDoParallel(cl)
                    cix <- foreach(i = 1:iterations, .combine=cbind) %dopar% {
                        tempMatrix <- apply(grx, 1, reg, PI = fit$PI[i, ])                           
                        tempMatrix           
                    }
                    stopCluster(cl)                
                    csx$x <- as.factor(csx$x)                
                    cat('DONE \n')
                    
                    for (i in 1:dim(cix)[2]) cix[,i]=quantile(ecdf(fit$Y[,2]),prob=exp(-1/cix[,i]))
                    csx$Zest <- apply(cix, 1, mean)
                    csx$lower <- apply(cix, 1, quantile, 0.025)
                    csx$upper <- apply(cix, 1, quantile, 0.975)
                    p <- p + 
                        geom_line(aes(y = Zest, group = as.factor(x)), size = 1.7) +
                        geom_ribbon(data=csx,aes(ymin = lower, ymax = upper), 
                                         linetype = 0, alpha = 0.1)
                }
                
                p <- p + geom_line(aes(y = Zest, group = as.factor(x)), size = 1.7)
                plot(p)
            }
            
        }
    } else if (which == 'qqplot') {
        k <- dim(fit$w)[1]
        iterations <- dim(fit$PI)[1]
        
        r <- array(dim = k)
        for (i in 1:k)
            r[i] <- qnorm(Gcond2D(fit$exceedances[i, 1], fit$exceedances[i, 2],
                                  fit$PIhat, 0))
        qq <- qqnorm(r, plot.it = FALSE)
        
        p <- ggplot(data = data.frame(x = qq$x, y = qq$y)) +
            geom_line(aes(x = qq$x, y = qq$y), lty = 2, col = "blue", size = 1.1) +
            geom_qq_line(aes(sample = r)) + 
            coord_fixed() +
            labs(y = "Sample quantiles", x = "Theoretical quantiles") +
            theme_classic() +
            theme(axis.text = element_text(size = 10),
                  axis.title = element_text(size = 10),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 10),
                  legend.key.size = unit(2, "line"))
        
        if (bands == TRUE) {
            cat('Learning...\n')            
            R <- matrix(nrow = k, ncol = iterations)
            traj <- matrix(nrow = iterations, ncol = k)
            for (i in 1:k)
                for (t in 1:iterations)
                    R[i, t] <- qnorm(Gcond2D(fit$exceedances[i, 1], fit$exceedances[i, 2],
                                             fit$PI[t, ], 0))
            for (t in 1:iterations)
                traj[t, ] <- qqnorm(R[, t], plot.it = FALSE)$y
            u.band <- apply(traj, 2, quantile, 0.975)
            l.band <- apply(traj, 2, quantile, 0.025)
            cat('DONE \n')
            
            p <- p +
                geom_ribbon(aes(x = qq$x, y = qq$y, ymin = l.band, ymax = u.band), 
                            linetype = 0, alpha = 0.1)
        }
        plot(p)
    } else if (which == 'qqboxplot') {
        k <- dim(fit$w)[1]
        iterations <- dim(fit$PI)[1]
        
        r <- array(dim = k)
        for (i in 1:k)
            r[i] <- qnorm(Gcond2D(fit$exceedances[i, 1], fit$exceedances[i, 2],
                                  fit$PIhat, 0))
        p <- ggplot(data = data.frame(y=r),aes(y=y)) + 
            qqboxplot::geom_qqboxplot(reference_dist="norm")
        plot(p)
    }
    
}

summary.regmani.BP <- function(fit, x = NULL, q = NULL) {
    ## Basic input validation
    if (sum(q >= 1) > 0  || sum(q < 0) > 0)
        stop('q must be inside the unit interval')
    
    cat('Learning...\n')  
    
    ## Regression manifold workhorse functions
    Gcond2D <- function(x, y, w, q) {
        ## conditional distribution estimate for Y | X = x  
        a1 <- w * fit$indices[, 1]
        a2 <- w * fit$indices[, 2]
        b1 <- 1 - pbeta(x / (x + y), fit$indices[, 1] + 1, fit$indices[, 2])
        b2 <- pbeta(x / (x + y), fit$indices[, 1], fit$indices[, 2] + 1)
        return(2 / fit$J * exp(-2 / fit$J * (x^{-1} * (a1%*%b1) + y^{-1} *
                                                 (a2%*%b2)) + x^{-1}) * (a1%*%b1) - q)
    }
    
    reg <- function(xq, PI) 
        uniroot(Gcond2D, c(0.0005, 10^12), x = xq[1],
                w = PI, q = xq[2], maxiter = 10000, check.conv = TRUE, tol = 1e-15)$root
    
    FX <- ecdf(fit$Y[, 1])
    gr_x <- expand.grid("x"=-1 / log(nrow(fit$Y) / (nrow(fit$Y) + 1) * FX(x)),"q"= q)
    
    cores <- detectCores()
    cat('using ', cores, ' cores...\n')            
    cl <- makeCluster(cores[1] - 1)
    print(cl)
    registerDoParallel(cl)
    
    ciRM <- foreach(i = 1:nrow(fit$PI), .combine=cbind) %dopar% {
        RM_ME <- apply(gr_x,1,reg,fit$PI[i,])
        temp <- quantile(ecdf(fit$Y[,2]),prob=exp(-1/RM_ME))
        temp
    }
    stopCluster(cl) 
    
    RM_pred <- rowMeans(ciRM)
    lq <- apply(ciRM,1,quantile,0.025)
    uq <- apply(ciRM,1,quantile,0.975)
    
    df <- data.frame(t(matrix(RM_pred,ncol=length(x))))
    df_lq <- data.frame(t(matrix(lq,ncol=length(x))))
    df_uq <- data.frame(t(matrix(uq,ncol=length(x))))
    rownames(df) <-  rownames(df_lq) <-  rownames(df_uq) <- apply(matrix(q*100,ncol=1),1,paste0,"%")
    colnames(df) <- colnames(df_lq) <- colnames(df_uq) <- apply(matrix(x,ncol=1),1,paste0)
    
    rtn <- list()
    rtn$table <- knitr::kable(df, escape=F, digits=4, 
                    caption = "Predicted quantiles of dependent variable (rows) \n evaluated for selected values of a covariate (columns)")
    rtn$lower <- knitr::kable(df_lq, escape=F, digits=4, 
                              caption = "lower bound of 95% credible interval")
    rtn$upper <- knitr::kable(df_uq, escape=F, digits=4, 
                              caption = "upper bound of 95% credible interval")
    return(rtn)
}
