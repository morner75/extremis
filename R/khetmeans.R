khetmeans <- function(y, centers, iter.max = 10, alpha = 0.5)
    UseMethod("khetmeans")

khetmeans.default <- function(y, centers, iter.max = 10, alpha = 0.5) {
    if(length(centers) > 1){k <- length(centers)}
    if(length(centers) == 1){k <- centers}
    N <- dim(y)[2] 
    T <- dim(y)[1]
    h <- 0.1 
    K <- floor((0.4258597) * T / (log(T))) 
    GRID <- seq(0, 1, length = 50) 
    
    Scf <- function(s, l) {
        y.ord <- sort(y[, l])
        I <- 1:T
        (1 / (K * h)) * sum((y[, l] > y.ord[T - K]) * G((s - (I / T)) / h))      
    }
    
    G <- function(x)
        ((15 / 16) * (1 - x^2)^2) * (x <= 1) * (-1 <= x)
    
    HILL <- numeric(N)
    for(l in 1: N) {
        y.ord <- sort(y[, l])
        HILL[l] <- (1 / K) *
            sum(log(y.ord[(T - K + 1):T]) - log(y.ord[T - K])) 
    }  
    
    ## Allocate cluster centers at random
    if(length(centers)>1){Ik <- centers}
    else{Ik <- sort(sample(seq(1, N, 1), k))}   
  
    mus <- function(s, j)
        Scf(s, Ik[j])      
    mugamma <- HILL[Ik]  
    M <- matrix(0, nrow = iter.max, ncol = N)
      
    ## sintegral is a function from the package Bolstad
    sintegral <- function (x, fx, n.pts = max(256, length(x))) 
    {
        if (class(fx) == "function") 
            fx = fx(x)
        n.x = length(x)
        if (n.x != length(fx)) 
            stop("Unequal input vector lengths")
        ap = approx(x, fx, n = 2 * n.pts + 1)
        h = diff(ap$x)[1]
        integral = h * (ap$y[2 * (1:n.pts) - 1] + 4 * ap$y[2 * (1:n.pts)] + 
                        ap$y[2 * (1:n.pts) + 1])/3
        results = list(value = sum(integral),
                       cdf = list(x = ap$x[2 * (1:n.pts)],
                                  y = cumsum(integral)))
        class(results) = "sintegral"
        return(results)
    }
    
    ## FIRST CYCLE
    cat(format("\nRendering\n"))
    cat(format("=========\n"))
    cat(1, "\n")
    clasification <- numeric(N)  
    for(i in 1:N) {    
        dist <- numeric(k)
        for(j in 1:k) {      
            V <- function(s) 
            (Scf(s, i) - mus(s, j))^2
            dist[j] <- alpha * (sintegral(GRID, Vectorize(V))$value) +
                (1 - alpha) * ((HILL[i] - mugamma[j])^2)
        }
        clasification[i] <- which.min(dist)
    }  
    mus.new <- function(s, j) {    
        aux.1 <- c()
        for(i in 1:N) 
            aux.1[i] <- Scf(s, i) * (1 * (clasification == j))[i]
        sum(aux.1) / sum(1 * (clasification == j))
    }
    
    mugamma.new <- numeric(k)  
    for(i in 1:k)
        mugamma.new[i] <- sum(HILL[which(clasification == i)]) /
            sum(1 * (clasification == i))  
    M[1, ] <- clasification
    
    ## CYCLES
    for(z in 2:iter.max) {
        cat(format(z), "\n")
        for(i in 1:N) {      
            dist <- numeric(k)      
            for(j in 1:k) {        
                V <- function(s)
                (Scf(s, i) - mus.new(s, j))^2
                dist[j] <- alpha * (sintegral(GRID, Vectorize(V))$value) +
                    (1 - alpha) * ((HILL[i] - mugamma.new[j])^2)
            }
            clasification[i] <- which.min(dist)
        }    
        mus.new <- function(s, j) {      
            aux.1 <- c()
            for(i in 1:N) 
                aux.1[i] <- Scf(s, i) * (1 * (clasification == j))[i]            
            sum(aux.1) / sum(1 * (clasification == j))
        }        
        mugamma.new <- numeric(k)    
        for(i in 1:k) 
            mugamma.new[i] <- sum(HILL[which(clasification == i)]) /
                sum(1 * (clasification == i))
        M[z, ] <- clasification    
        if(prod(1 * (M[z, ] == M[(z - 1), ])) == 1)
            break          
    }
    
    ## Organize and return outputs    
    outputs <- list(mugamma.new = mugamma.new, Y = y, n.clust = k,
                    mus.new = mus.new, clusters = M[z - 1, ])
    class(outputs) <- "khetmeans"
    return(outputs)
}


plot.khetmeans <- function(x, c.c = FALSE, xlab = "w",
                           ylab = "Scedasis function", ...) {  
    N <- dim(x$Y)[2] 
    T <- dim(x$Y)[1]
    h <- 0.1 
    K <- floor((0.4258597) * T / (log(T))) 
    grid <- seq(0, 1, length = 100)
    
    Scf <- function(s, l) {
        y.ord <- sort(x$Y[, l])
        I <- 1:T
        (1 / (K * h)) * sum((x$Y[, l] > y.ord[T - K]) * G((s - (I / T)) / h))      
    }
    
    G <- function(x)
        ((15 / 16) * (1 - x^2)^2) * (x <= 1) * (-1 <= x)
    
    s1 <- lapply(grid, Scf, 1)
    plot(grid, s1, type = 'l', xlab = xlab, ylab = ylab, col = 'gray', ...)
    
    for(i in 1:N) {
        s <- lapply(grid, Scf, i)
        lines(grid, s, type = 'l', lwd = 8, col = 'gray')
    }     
    
    if(c.c) {
        for(j in 1:x$n.clust) {
            y <- lapply(grid, x$mus.new, j)
            lines(grid, y, type = 'l', lwd = 8, col = "blue", lty = 6)
        }
    } 
}
