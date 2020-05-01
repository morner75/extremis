\encoding{UTF-8}
\name{angdensity}
\alias{angdensity}
\alias{angdensity.default}
\title{Empirical-Likelihood Based Inference for the Angular Density}
\description{
  This function computes empirical-likelihood based estimators for the
  angular distribution function of a bivariate extreme value
  distribution.
}
\usage{
angdensity(Y, tau = 0.95, nu, grid = seq(0.01, 0.99, length = 2^8),
	   method = "euclidean", raw = TRUE)
}
\arguments{
  \item{Y}{data frame with two columns from which the estimate is to
  be computed.}
  \item{tau}{value used to threshold the data; by default it is set as
  the 0.95 quantile of the pseudo-radius.}
  \item{nu}{concentration parameter of beta distribution which
  controls the amount of smoothing.}
  \item{grid}{grid with coordinates of the points where the angular
  density is estimated; by default \code{grid = seq(0.01, 0.99, length
  = 2^8)}.}
  \item{method}{a character string setting the method to be used. By
  default \code{method = "euclidean"}, the other option being
  \code{method = "empirical"}. See details.}
  \item{raw}{logical; if \code{TRUE}, \code{Y} will be converted to
  unit Fréchet scale. If \code{FALSE}, \code{Y} will be understood as
  already in unit Fréchet scale.}
}
\value{
  \item{h}{the estimate angular density values.}
  \item{grid}{grid with coordinates of the points where the angular
  density is estimated.}
  \item{w}{pseudo-angles.}
  \item{nu}{concentration parameter of the Beta-kernel.}
  \item{Y}{raw data.}
  The \code{plot} method depicts the smooth angular density.
}
\details{
  The smooth angular density was introduced in by de Carvalho et al
  (2013). \code{method = "euclidean"} implements the version of the
  method based on Euclidean likelihood weights, whereas \code{method =
  "empirical"} uses Empirical likelihood weights. 
}
\references{
  de Carvalho, M., Oumow, B., Segers, J. and Warchol, M. (2013) A
  Euclidean likelihood estimator for bivariate tail dependence.
  \emph{Communications in Statistics---Theory and Methods}, 42,
  1176--1192.
}
\examples{
## de Carvalho et al (2013, Fig. 7)
data(beatenberg)
attach(beatenberg)
fit <- angdensity(beatenberg, tau = 0.98, nu = 163, raw = FALSE)
plot(fit)
rug(fit$w)
}
\author{Miguel de Carvalho}
