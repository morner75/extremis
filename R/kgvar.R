kgvar <- function(y, centers, iter.max = 10, conf.level = 0.95)
    UseMethod("kgvar")

kgvar.default <- function(y, centers, iter.max = 10, conf.level = 0.95) {    
    if(length(centers) > 1){k <- length(centers)}
    if(length(centers) == 1){k <- centers}
    N <- dim(y)[2] 
    T <- dim(y)[1]
    h <- 0.1 
    p <- conf.level
    K <- floor((0.4258597) * T / (log(T))) 
    GRID <- seq(0.05, 1, length = 50) 
    I <- 1:T
    
    HILL <- numeric(N)
    
    G <- function(x)
        ((15 / 16) * (1 - x^2)^2) * (x <= 1) * (-1 <= x)
    
    a <- numeric(N)  
    for(l in 1:N) {
        y.ord <- sort(y[, l])
        U <- y.ord[T - K]
        HILL[l] <- (1 / K) * sum(log(y.ord[(T - K + 1):T]) - log(y.ord[T - K])) 
       
        SCEDASIS <- function(s)
        (1 / (K * h)) * sum((y[, l] > y.ord[T - K]) * G((s - (I / T)) / h))    
        
        SS <- numeric(T)    
        for(j in 1:T)
            SS[j] <- (SCEDASIS(j / T))^(1 / HILL[l])  
        a[l] <- (U * (K)^(HILL[l])) / (sum(SS))^(HILL[l])    
    }
    
    Varf <- function(s, l) {
        y.ord <- sort(y[, l])
        log((a[l]/((1-p)^(HILL[l])))*(1 / (K * h)) *
            sum((y[, l] > y.ord[T - K]) * G((s - (I / T)) / h)))     
    }
    
    if(length(centers) > 1){Ik = centers}  
    if(length(centers) == 1){Ik <- sort(sample(seq(1, N, 1), k))}   
    
    mus <- function(s, j)
        Varf(s, Ik[j])      
  
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
            V<-function(x) 
            (Varf(x, i) - mus(x, j))^2      
            dist[j] <- sintegral(GRID, Vectorize(V))$value      
        }
        clasification[i] <- which.min(dist)
    }  
    
    mus.new <- function(s, j) {    
        aux.1 <- c()
        for(i in 1:N) 
            aux.1[i] <- Varf(s, i) * (1 * (clasification == j))[i]
        sum(aux.1) / sum(1 * (clasification == j))
    }
    
    M[1, ] <- clasification
    
    ## CYCLES
    for(z in 2:iter.max) {
        cat(format(z), "\n")
        for(i in 1:N) {      
            dist <- numeric(k)      
            for(j in 1:k) {        
                V.new<-function(s)
                (Varf(s, i) - mus.new(s, j))^2
                dist[j] <- sintegral(GRID, Vectorize(V.new))$value        
            }
            clasification[i] <- which.min(dist)
        }    
        mus.new <- function(s, j) {      
            aux.1 <- c()
            for(i in 1:N) 
                aux.1[i] <- Varf(s, i) * (1 * (clasification == j))[i]            
            sum(aux.1) / sum(1 * (clasification == j))
        }        
        
        M[z, ] <- clasification    
        if(prod(1 * (M[z, ] == M[(z - 1), ])) == 1)
            break          
    }
    
    ## Organize and return outputs    
    outputs <- list(Y = y, n.clust = k, scale.param = a,
                    conf.level = conf.level, hill = HILL, var.new = mus.new,
                    clusters = M[z - 1,])
    class(outputs) <- "kgvar"
    return(outputs)
}


plot.kgvar <- function(x, c.c = FALSE, xlab = "w",
                       ylab = "Value-at-risk function", ...) {  
    N <- dim(x$Y)[2] 
    T <- dim(x$Y)[1]
    h <- 0.1 
    a <- x$scale.param
    p <- x$conf.level
    HILL <- x$hill
    K <- floor((0.4258597) * T / (log(T))) 
    grid <- seq(0, 1, length = 100)
    
    Varf <- function(s, l) {
        y.ord <- sort(x$Y[, l])
        I <- 1:T
        (a[l]/((1-p)^(HILL[l])))*(1 / (K * h)) *
            sum((x$Y[, l] > y.ord[T - K]) * G((s - (I / T)) / h))     
    }
    
    G <- function(x)
        ((15 / 16) * (1 - x^2)^2) * (x <= 1) * (-1 <= x)
  
    s1 <- lapply(grid, Varf, 1)
    plot(grid, s1, type = 'l', xlab = xlab, ylab = ylab, col = 'gray', ...)
  
    for(i in 1:N) {
        s <- as.numeric(lapply(grid,Varf,i))
        lines(grid, s, type = 'l', lwd = 8, col = 'gray')
    }     
  
    Newvar2 <- function(s,j){exp(x$var.new(s,j))}    
    if(c.c) {    
        for(j in 1:x$n.clust){
            y <- lapply(grid, Newvar2, j)
            lines(grid, y, type = 'l', lwd = 8, col = "blue", lty = 6)
        }   
    }  
}
