\encoding{UTF-8}
\name{bevfd.kernel}
\alias{bevfd.kernel}
%% \alias{bev.kernel.default}
\title{Kernel Inference for Bivariate Extremes Conditional on
  Functional Covariates}
\description{
  This function computes kernel-based estimators of the conditional
  angular density and angular distribution function of a bivariate
  extreme value distribution conditionally on a functional covariate.
}
\usage{
bevfd.kernel(Y, Xt, X, tau = 0.95, grid = seq(0.01, 0.99, length = 2^8),
   smoothing = "lscv", raw = TRUE, target = "pdf")
}
\arguments{
  \item{Y}{data frame with two columns from which the estimate is to
    be computed.}
  \item{Xt}{matrix of size \code{n x m} with discrete evaluation of
    functional covariates, where \code{n} is the number of curves and
    \code{m} is the number of points per curve.}
  \item{X}{function conditionally on which the angular measure estimate
    is evaluated.}
   \item{tau}{value used to threshold the data; by default it is set as
     the 0.95 quantile of the pseudo-radius.}
   \item{smoothing}{by default, \code{smoothing = 'lscv'}; otherwise, a
     list containing the smoothing parameters \code{b}, \code{nu}, and
     \code{c} can be set.}
   \item{grid}{grid with coordinates of the points where the angular
     density is estimated; by default \code{grid = seq(0.01, 0.99, length
       = 2^8)}.}
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
  \item{grid}{grid with coordinates of the points where the angular
  density is estimated.}
  \item{w}{pseudo-angles.}
  \item{smoothing}{list with smoothing parameters.}
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
Alabdulathem, A.,de Carvalho, M., and Leonelli, M. (2022). Functional Covariate-Adjusted Extremal Dependence.
}
\examples{
## ================= ## 
##  Simulated data   ##
## ================= ## 
set.seed(1234)
n <- 1000 
Xt <- Sim.DiffProc::BM(M = n)
Xt <- fda.usc::fdata(as.data.frame(Xt))
alpha <- cos(fda.usc::norm.fdata(Xt))^2
for (i in 1:n) Y <- evd::rmvevd(n, dep = alpha[i], model = "log", d = 2)
## set curve conditionally on which the fit is plotted
t <- seq(0.01, 0.99, length.out = n)
X <- fda.usc::fdata(sin(3 * t))
## fit conditional angular measure
fit <- bevfd.kernel(Y = Y, X = X, Xt = Xt,
                    smoothing = list(b = 2, nu = 1, c = 15))
## to plot against the true
alpha1 <- cos(fda.usc::norm.fdata(X))^2
###############################
require(ggplot2)
plot(fit)+ geom_line( aes( y = evd::hbvevd(fit$grid, dep =  alpha1[1,], 
                model = "log", half = TRUE)))

}
