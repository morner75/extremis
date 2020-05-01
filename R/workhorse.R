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

##  ========================================================================  ##
##   'Workhorse' function for plotting pointwisely approximated trajectories  ##
##  ------------------------------------------------------------------------  ##
plot.workhorse <- function(grid, traj, trajhat, square, bands = TRUE,
                           last = FALSE, xxlab, yylab, xlim = NULL,
                           ylim = NULL, ...) {
    xlim <- if (is.null(xlim)) 
        range(grid)
    else xlim

    xxlab <- if (is.null(xxlab))
        NULL 
    else xxlab

    if (square == TRUE)
        par(pty = "s")
    
    if (last == TRUE) {
        T <- nrow(traj)
        if (T < 100)
            stop('"last = TRUE" requires more than 100 trajectories')
        ylim <- if (is.null(ylim)) 
            c(0, max(traj[(T - 100 + 1):T, ]))
        else ylim
        
        plot(grid[1], trajhat[1], xlim = xlim, ylim = ylim, type = "l",
             xlab = xxlab, ylab = yylab, ...)
        
        lines(grid, trajhat, type = "l", col = "black", lwd = 1, ...)
        for (t in 1:100)
            lines(c(0, grid), c(0, traj[(T - t + 1), ]), type = "l",
                  lwd = 1, col = t)
        bands = FALSE
    }

    if (bands == TRUE) {
        l.band <- apply(traj, 2, quantile, probs = 0.025)
        u.band <- apply(traj, 2, quantile, probs = 0.975)
        
        ylim <- if (is.null(ylim)) 
                    c(0, max(u.band))
        else ylim

        plot(grid[1], trajhat[1], xlim = xlim, ylim = ylim, type = "l",
             xlab = xxlab, ylab = yylab, ...)
        polygon(x = c(rev(grid),grid), y = c(rev(l.band), u.band),
                border = NA, col = "lightgray")
        lines(c(0, grid), c(0, trajhat), type = "l", col = "black",
              lwd = 3, ...)
    }
    
    if (bands == FALSE & last == FALSE)
        plot(grid, trajhat, xlim = xlim, ylim = ylim, type = "l", lwd = 3,
             xlab = xxlab, ylab = yylab, ...)
}
