\encoding{UTF-8}
\name{khetmeans}
\alias{khetmeans}
\alias{khetmeans.default}
\title{K-Means Clustering for Heteroscedastic Extremes}
\description{
  This function performs k-means clustering for heteroscedastic extremes.
}
\usage{
  khetmeans(y, centers, iter.max = 10, alpha = 0.5)
}
\arguments{
  \item{y}{data frame from which the estimate is to be computed; first
    column corresponds to time and the second to the remainder of interest.} 
  \item{centers}{the number of clusters or a set of initial (distinct)
    cluster centres. If a number, a random set of (distinct) rows in
    \code{y} is chosen as the initial centers.}
  \item{iter.max}{the maximum number of iterations allowed. The
   default is 10.}
  \item{alpha}{the tuning parameter. The default is 0.5.}
}
\value{khetmeans returns an object of class "\code{khetmeans}" which has a
  fitted method. It is a list with at least the following components:
  \item{mus.new}{cluster center scedasis density.}
  \item{mugamma.new}{cluster center extreme value index.}
  \item{clusters}{cluster allocation.}
  \item{Y}{raw data.}
  \item{n.clust}{number of clusters.}
  The \code{plot} method depicts the k-means clustering for
  heteroscedastic extremes. If \code{c.c} is \code{TRUE}, the method
  displays the cluster means.
}
\details{The intermediate sequence \eqn{\kappa_T} is chosen
  proportional to \eqn{T/\log T}.
}
\references{Rubio, R., de Carvalho, M. and Huser, R. (2018) 
  Similarity-Based Clustering of Extreme Losses from the London Stock
  Exchange. Submitted.}
\examples{
\dontrun{
## Example 1  (Scenario B, T = 5000)
## This example requires package evd 
require(evd)
set.seed(12)
T <- 5000
n <- 30
b <- 0.1
gamma1 <- 0.7
gamma2 <- 1
grid <- seq(0, 1, length = 100)
c2 <- function(s)
    dbeta(s, 2, 5)
c3 <- function(s)
    dbeta(s, 5, 2)
X <- matrix(0, ncol = T, nrow = n)
for(i in 1:5)
  for(j in 1:T)
    X[i,  j] <- rgev(1, c2(j / T), c2(j / T), gamma1)
for(i in 6:15)
  for(j in 1:T)
    X[i,  j] <- rgev(1, c2(j / T), c2(j / T), gamma2)
for(i in 16:20)
  for(j in 1:T)
    X[i,  j] <- rgev(1, c3(j / T), c3(j / T), gamma1)
for(i in 21:30)
  for(j in 1:T)
    X[i,  j] <- rgev(1, c3(j / T), c3(j / T), gamma2)
Y <- t(X)
fit <- khetmeans(Y, centers = 4)
plot(fit, c.c = TRUE)
lines(grid, c2(grid), type = 'l', lwd = 8, col = 'black')
lines(grid, c3(grid), type = 'l', lwd = 8, col = 'black')
}

\dontrun{
## Example 2 (Overlapping version of Fig. 5 in Supplementary Materials)
data(lse)
attach(lse)
y <- -apply(log(lse[, -1]), 2, diff)
fit <- khetmeans(y, centers = 3)
plot(fit, c.c = TRUE, ylim = c(0, 3))
}
}
\author{Miguel de Carvalho, Rodrigo Rubio.}