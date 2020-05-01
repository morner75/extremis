\encoding{UTF-8}
\name{angcdf}
\alias{angcdf}
\alias{angcdf.default}
\title{Empirical-Likelihood Based Inference for the Angular Measure}
\description{
  This function computes empirical-likelihood based estimators for the
  angular distribution function of a bivariate extreme value
  distribution.
}
\usage{
angcdf(Y, tau = 0.95, method = "euclidean", raw = TRUE)
}
\arguments{
  \item{Y}{data frame with two columns from which the estimate is to
  be computed.}
  \item{tau}{value used to threshold the data; by default it is set as
  the 0.95 quantile of the pseudo-radius \code{tau = 0.95}.}
  \item{method}{a character string setting the method to be used. By
  default \code{method = "euclidean"}, the other option being
  \code{method = "empirical"}. See details.}
  \item{raw}{logical; if \code{TRUE}, \code{Y} will be converted to
  unit Fréchet scale. If \code{FALSE}, \code{Y} will be understood as
  already in unit Fréchet scale.}
}
\value{
  \item{H}{angular distribution function.}
  \item{w}{pseudo-angles.}
  \item{Y}{data.}
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
fit <- angcdf(beatenberg, tau = 0.98, raw = FALSE)
plot(fit)
rug(fit$w)
}
\author{Miguel de Carvalho}
