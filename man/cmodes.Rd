\encoding{UTF-8}
\name{cmodes}
\alias{cmodes}
\alias{cmodes.default}
\title{Mode Mass Function}
\description{
  This function computes the mode mass function.
}
\usage{
cmodes(Y, thresholds = apply(Y[, -1], 2, quantile, probs =
                 0.95), nu = 100, ...)
}

\arguments{
  \item{Y}{data frame from which the estimate is to be computed; first
    column corresponds to time and the second to the variable of interest.}  
  \item{thresholds}{values used to threshold the data \code{y}; by default
    \code{threshold = quantile(y, 0.95)}.}
  \item{nu}{concentration parameter of beta kernel used to smooth mode
  mass function.}
  \item{...}{further arguments for \code{density} methods.}
}

\value{
  \item{c}{scedasis density estimators.}
  \item{k}{number of exceedances above the threshold.}
  \item{w}{standardized indices of exceedances.}
  \item{Y}{raw data.}
  The \code{plot} method depicts the smooth mode mass function along
  with the smooth scedasis densities. 
}
\details{
  The scedasis functions on which the mode mass function is based are
  computed using the default \code{"nrd0"} option for bandwidth.
}
\references{
  Rubio, R., de Carvalho, M., and Huser, R. (2018)
  Similarity-Based Clustering of Extreme Losses from the London Stock
  Exchange. Submitted. 
}
\examples{
data(lse)
attach(lse)
nlr <- -apply(log(lse[, -1]), 2, diff)
Y <- data.frame(DATE[-1], nlr)
T <- dim(Y)[1]
k <- floor((0.4258597) * T / (log(T)))
fit <- cmodes(Y, thresholds = as.numeric(apply(nlr, 2, sort)[T - k, ]),  
              kernel = "biweight", bw = 0.1 / sqrt(7), nu = 100)
plot(fit)
}
\author{Miguel de Carvalho}