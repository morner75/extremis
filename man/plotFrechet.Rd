\encoding{UTF-8}
\name{plotFrechet}
\alias{plotFrechet}
\alias{plotFrechet.default}
\title{Unit Fréchet Scatterplot in Log-log Scale}
\description{
  This function depicts a scatterplot of bivariate data transformed 
  to unit Fréchet scale. 
}
\usage{
plotFrechet(Y, tau = 0.95, raw = TRUE, ...)
}
\arguments{
  \item{Y}{list with data from which the estimates are to be computed.}
  \item{tau}{value used to threshold the data \code{y}; by default
    \code{treshold = quantile(y, 0.95)}.}
  \item{raw}{logical; if \code{TRUE}, \code{Y} will be converted to unit Fréchet
  scale. If \code{FALSE}, \code{Y} will be understood as already in unit Fréchet scale.}
  \item{...}{other arguments to be passed to \code{plot}.}
}
\details{
  The solid line corresponds to the boundary threshold in the log-log
  scale, with both axes being logarithmic.
}
\examples{
## de Carvalho et al (2013, Fig. 5)
data(beatenberg)
plotFrechet(beatenberg, xlab = "Forest Cover", ylab = "Open Field",
            raw = FALSE)
}
\author{Miguel de Carvalho}
