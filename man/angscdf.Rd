\encoding{UTF-8}
\name{angscdf}
\alias{angscdf}
\alias{angscdf.default}
\title{Smooth Empirical-Likelihood Based Inference for the Angular Measure}
\description{
  This function computes smooth empirical-likelihood based estimators for the
  angular distribution function of a bivariate extreme value
  distribution.
}
\usage{
angscdf(Y, tau = 0.95, nu, grid = seq(0.01, 0.99, length = 2^8),
	method = "euclidean", raw = TRUE)
}
\arguments{
  \item{Y}{data frame with two columns from which the estimate is to
  be computed.}
  \item{tau}{value used to threshold the data; by default it is set as
  the 0.95 quantile of the pseudo-radius \code{tau = 0.95}.}
  \item{nu}{concentration parameter of beta distribution which
  controls the amount of smoothing.}
  \item{grid}{grid with coordinates of the points where the angular
  measure is estimated; by default \code{grid = seq(0.01, 0.99, length
  = 2^8)}.}
  \item{method}{a character string setting the method to be used. By
  default \code{method = "euclidean"}, the other option being
  \code{method = "empirical"}. See details.}
  \item{raw}{logical; if \code{TRUE}, \code{Y} will be converted to
  unit Fréchet scale. If \code{FALSE}, \code{Y} will be understood as
  already in unit Fréchet scale.}
}
\value{
  \item{H}{the estimated angular distribution function values.}
  \item{grid}{grid with coordinates of the points where the angular
  measure is estimated.}
  \item{w}{pseudo-angles.}
  \item{nu}{concentration parameter of the Beta-kernel.}
  \item{Y}{raw data.}
  The \code{plot} method depicts the empirical likelihood-based
  angular distribution function.
}
\details{
  \code{method = "euclidean"} implements the maximum Euclidean
  likelihood spectral distribution function as introduced by de
  Carvalho et al (2013). \code{method = "empirical"} implements the
  maximum Empirical likelihood spectral distribution function as
  introduced by Einmahl and Segers (2009). 
}
\references{
  de Carvalho, M., Oumow, B., Segers, J. and Warchol, M. (2013) A
  Euclidean likelihood estimator for bivariate tail dependence.
  \emph{Communications in Statistics---Theory and Methods}, 42,
  1176--1192.

  Einmahl, J. H. J., and Segers, J. (2009) Maximum empirical
  likelihood estimation of the spectral measure of an extreme-value
  distribution.  \emph{The Annals of Statistics}, 37, 2953--2989.
}
\examples{
## de Carvalho et al (2013, Fig. 7)
data(beatenberg)
attach(beatenberg)
fit <- angscdf(beatenberg, tau = 0.98, nu = 163, raw = FALSE)
plot(fit)
rug(fit$w)
}
\author{Miguel de Carvalho}