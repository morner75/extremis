\encoding{UTF-8}
\name{bev.kernel}
\alias{bev.kernel}
%% \alias{bev.kernel.default}
\title{Kernel Empirical-Likelihood Inference for Bivariate Extremes}
\description{
  This function computes kernel empirical-likelihood based estimators of
  the angular density and angular distribution function of a bivariate 
  extreme value distribution.
}
\usage{
bev.kernel(Y, tau = 0.95, nu, grid = seq(0.01, 0.99, length = 2^8),
           method = "euclidean", raw = TRUE, target = "pdf")
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
  \code{method = "empirical"}; see details.}
  \item{raw}{logical; if \code{TRUE}, \code{Y} will be converted to
  unit Fréchet scale. If \code{FALSE}, \code{Y} will be understood as
  already in unit Fréchet scale.}
  \item{target}{value set to \code{"pdf"}, \code{"cdf"}, or
    \code{"Pickands"} depending on whether the goal is to estimate an
    angular density, angular distribution function, or Pickands
    dependence function, respectively; by default: \code{target =
      "pdf"}.}
}
\value{
  \item{estimate}{estimated target of interest evaluated over
    \code{grid}.}
  \item{p}{Euclidean likelihood weights, if \code{method = "euclidean"},
    or empirical likelihood weights, if \code{method = "emprical"}.}
  \item{grid}{grid with coordinates of the points where the angular
  density is estimated.}
  \item{w}{pseudo-angles.}
  \item{nu}{concentration parameter of the Beta kernel.}
  \item{target}{string with target of interest.}
  The \code{plot} method depicts the estimated object [angular density
  (\code{target = 'pdf'}), angular distribution function (\code{target =
    'cdf'}) or Pickands dependence function (\code{target =
    'Pickands'})].
}
\details{
  The smooth angular density was introduced by de Carvalho et al
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
## ================= ## 
##  Beatenberg data  ##
## ================= ## 
## de Carvalho et al (2013, Fig. 7)
data(beatenberg)
fit <- bev.kernel(beatenberg, tau = 0.98, nu = 163, raw = FALSE)
plot(fit)

## To confirm that the produced estimate obeys the normalization and moment
## constraints, type: 
sum(fit$p)            # = 1 
sum(fit$w * fit$p)    # = 1/2 

\dontrun{
## Run the lines below to compare with the Bernstein polynomial angular density
## estimator of Hanson et al (2017)
fit.BP <- bev.BP(beatenberg, tau = 0.98, raw = FALSE)
chart <- seq(0.01, 0.99, by = 0.01)
lines(fit$grid, fit$trajhat)

## Run the lines below to compare with the Bernstein polynomial angular measure
## estimator of de Carvalho et al (2017)
fit.BP <- bev.BP(beatenberg, tau = 0.98, raw = FALSE, target = "cdf")
plot(fit.BP)
rug(fit$w)
fit.kernel <- bev.kernel(beatenberg, tau = 0.98, nu = 163,
                         raw = FALSE, target = "cdf")
lines(fit.kernel$grid, fit.kernel$estimate)

## de Carvalho et al (2013, Fig. 8)
fit <- bev.kernel(beatenberg, tau = 0.98, nu = 163, raw = FALSE,
                  target = "Pickands")
plot(fit)
}
}
