\encoding{UTF-8}
\name{cdf}
\alias{cdf}
\alias{cdf.default}
\title{Empirical Scedasis Distribution Function}
\description{
  This function computes the empirical scedasis distribution function.
}
\usage{
cdf(Y, threshold = quantile(Y[, 2], 0.95))
}

\arguments{
  \item{Y}{data frame from which the estimate is to be computed; first
  column corresponds to time and the second to the variable of interest.}
  \item{threshold}{value used to threshold the data \code{y}; by default
    \code{threshold = quantile(Y[, 2], 0.95)}.}
}

\value{
  \item{C}{empirical scedasis distribution function.}
  \item{w}{standardized indices of exceedances.}
  \item{k}{number of exceedances above a threshold.}  
  \item{Y}{raw data.}
  The \code{plot} method depicts the empirical cumulative scedasis
  function, and the reference line for the case of constant frequency of
  extremes over time (if \code{uniform = TRUE}). 
}
\details{
  The empirical scedasis distribution function was introduced by Einmahl
  et al (2016).
}

\references{
  Einmahl, J. H., Haan, L., and Zhou, C. (2016) Statistics of
  heteroscedastic extremes. \emph{Journal of the Royal Statistical
  Society: Ser. B}, 78(1), 31--51.
}
\examples{
data(sp500)
attach(sp500)
Y <- data.frame(date[-1], -diff(log(close)))
fit <- cdf(Y)
plot(fit)
plot(fit, original = FALSE)
}
\author{Miguel de Carvalho}
