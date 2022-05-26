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

regmani <- function(model = "dependence", dep, asy = c(1, 1), alpha, beta,
                    q = NULL, pred = NULL,...)
    UseMethod("regmani")

regmani.default <- function(model = "dependence", dep, asy = c(1, 1), alpha,
                            beta, q = NULL, pred = NULL, ...) {    
    ## Basic input validation
    if (sum(model == c("dependence", "independence", "log", "alog", "hr",
                       "neglog", "aneglog", "bilog", "negbilog", "ct",
                       "amix")) != 1)
        stop('model must be one of "dependence", "independence", "log", "alog",
             "hr", "neglog", "aneglog", "bilog", "negbilog", "ct" or "amix"')

    xmin <- 0.005
    xmax <- 20
    qmin <- 0.01
    qmax <- 0.99
    
    if (is.null(q) == TRUE & is.null(pred) == TRUE) {
        target <- 'full'
        l <- 30 
        grid <- expand.grid(x = seq(xmin, xmax, length = l),
                            q = seq(0.01, 0.99, length = l))            
    } else if (is.numeric(q) == TRUE) {
        ## input validation
        if (sum(q >= 1) > 0  || sum(q < 0) > 0)
            stop('q must be inside the unit interval')
        target <- 'fixed_q'
        l <- 50
        grid <- expand.grid(x = seq(xmin, xmax, length = l), q = q)
    } else if (is.numeric(pred) == TRUE) {
        ## input validation
        if (sum(pred <= 0) > 0) 
            stop('"pred" must be positive')
        target <- 'fixed_x'
        l <- 50
        grid <- expand.grid(x = pred, q = seq(qmin, qmax, length = l))
    }
    
    if (model == "dependence") {
        rm <- data.frame(x = grid[, 1], q = grid[, 2], z = grid[, 1])
    } else if (model == "independence") {
        rm <- data.frame(x = grid[, 1], q = grid[, 2], z = -1 / log(grid[, 2]))
    } else if (sum(model == c("log", "hr", "neglog")) == 1 & missing(asy)) {
        reg <- function(xq, model, dep, alpha, beta) {
            Gcond2D <- function(y, xq, model, dep, alpha, beta) 
                evd::ccbvevd(matrix(c(exp(-1 / y), exp(-1 / xq[1])), nrow = 1,
                                    ncol = 2), mar = 2, model = model, dep = dep,
                             alpha = alpha, beta = beta) - xq[2]        
            arg <- uniroot(Gcond2D, c(0.005, 10^10), xq = xq, model = model,
                           dep = dep, alpha = alpha, beta = beta)$root
            return(arg)
        }
        rm <- data.frame(x = grid[, 1], q = grid[, 2],
                         z = apply(grid, 1, reg, model, dep, alpha, beta))
        
    } else if (model == "ct") {
        reg <- function(xq, model, alpha, beta) {
            Gcond2D <- function(y, xq, model, alpha, beta) {
                qq = (alpha / y) / (alpha / y + beta / xq[1])
                gamma = (alpha + beta + 2) * (alpha + beta + 1)
                evd::pbvevd(c(xq[1],y), model = model, mar1 = c(1, 1, 1),
                            mar2 = c(1, 1, 1), alpha = alpha, beta = beta) * 
                    exp(1 / xq[1]) * (pbeta(qq, alpha + 1, beta, lower.tail = F) +
                                      (alpha + 1) * beta / gamma * dbeta(qq, alpha + 2, beta + 1) -
                                      xq[1] / y * alpha * (beta + 1) / gamma *
                                      dbeta(qq, alpha + 1, beta + 2)) - xq[2]
            } 
            arg <- uniroot(Gcond2D, c(0.005, 10^10), xq = xq, model = model,
                           alpha = alpha, beta = beta)$root
            return(arg)
        }
        rm <- data.frame(x = grid[, 1], q = grid[, 2],
                         z = apply(grid, 1, reg, model, alpha, beta))        
    } else if (sum(model == c("alog", "aneglog")) == 1 ){
        reg <- function(xq, model, dep, alpha, beta, asy) {
            Gcond2D <- function(y, xq, model, dep, alpha, beta, asy) 
                evd::ccbvevd(matrix(c(exp(-1 / y), exp(-1 / xq[1])), nrow = 1,
                                    ncol = 2), mar = 2, model = model, dep = dep,
                             alpha = alpha, beta = beta) - xq[2]        
            arg <- uniroot(Gcond2D, c(0.005, 10^10), xq = xq, model = model,
                           dep = dep, alpha = alpha, beta = beta, asy = asy)$root
            return(arg)
        }
        rm <- data.frame(x = grid[, 1], q = grid[, 2],
                         z = apply(grid, 1, reg, model, dep, alpha, beta, asy = asy))
    }       
    
    ## Organize and return outputs
    outputs <- list(rm = rm, target = target)
    class(outputs) <- "regmani"
    return(outputs)
}

plot.regmani <- function(chart, cex = 1, ...) {
    if (chart$target == 'full') {
        cex <- 1
        ## Setting transparency to every color in the vector
        blues_pal <- colorRampPalette(c("darkblue", 'blue', 'cyan',
                                        'green4', 'yellow'), bias = 1.5)(500)    
        blues_pal_transp <- transparent(blues_pal, trans.val = 0.4)
        
        wireframe(z ~ x * q, data = chart$rm,
                  scales = list(arrows = FALSE, col = "black", cex = cex),
                  xlab = list("x", cex = cex),
                  ylab = list("q", cex = cex),
                  zlab = list(latex2exp::TeX('$y_{q | x}$'), rot = 90, cex = cex),
                  cex = cex,
                  par.settings = list(axis.line = list(col = "transparent")),
                  drape = T, colorkey = F, par.box = list(col = 1, lty = 1),
                  col.regions = blues_pal_transp,
                  screen = list(z = 120, x = -80, y = 5), ...)        
    } else if (chart$target == 'fixed_q') {
        csp <- data.frame(x = chart$rm$x, q = chart$rm$q, z = chart$rm$z)
        ggplot(data = csp, aes(x = x, y = z, col = as.factor(q))) +
            geom_line(aes(x = x, y = z, group = as.factor(q)), size = 1.1) +
            theme_classic() +
            coord_fixed(xlim = c(0, 20), ylim = c(0, 20)) +
            labs(y = "y", x = "x") +
            guides(col = guide_legend(title = "q", order = 2), 
                   linetype = guide_legend(title = "Regression line", order = 1)) + 
            theme(axis.text = element_text(size = 10), 
                  axis.title = element_text(size = 10),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 10), 
                  legend.key.size = unit(2, "line"))        
    } else if (chart$target == 'fixed_x') {
        csx <- data.frame(x = chart$rm$x, q = chart$rm$q, z = chart$rm$z)
        ggplot(data = csx, aes(x = q, y = z, col = as.factor(x))) +
            geom_line(aes(x = q, y = z, group = as.factor(x)), size = 1.1) +
            theme_classic() +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, max(csx$z, csx$z))) +
            labs(y = "Conditional quantile", x="q") +
            guides(linetype = guide_legend(title = " "),
                   col = guide_legend(title = "x")) + 
            theme(axis.text=element_text(size = 10),
                  axis.title = element_text(size = 10),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 10),
                  legend.key.size = unit(2,"line"))
        
    }
}
