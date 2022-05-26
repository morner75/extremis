\encoding{UTF-8}
\name{chit}
\alias{chit}
\alias{chibart}
\title{Bayesian P-Spline Modeling of Time-Varying Extremal Dependence}
\description{
  These functions compute time-varying extremal dependence measures via
  Bayesian P-splines. 
}
\usage{
chit(XY, tau = 0.95, prior = list(a = 0.001, b = 0.001),
     T = 10000, burn = 5000, link = 'probit', knot = 20, dif = 2, 
     raw = TRUE)
chibart(XY, tau = 0.95, prior = list(a = 0.001, b = 0.001),
     T = 10000, burn = 5000, link = 'probit', knot = 20, dif = 2, 
     raw = TRUE)
}
\arguments{
  \item{XY}{data frame from which the estimate is to be computed; first
    column corresponds to time, the second and third columns correspond
    to the variables of interest \code{X} and \code{Y} respectively.}
  \item{tau}{value used to threshold the data \code{z = min(X, Y)}; by default
    \code{threshold = quantile(z, tau)}.}
  \item{prior}{list with hyperparameters for inverse Gamma prior for the
    variances.}
  \item{T}{number of MCMC iterations; by default \code{T = 10000}.}
  \item{burn}{number of burn in iterations; by default \code{burn = 5000}.}
  \item{link}{a character string giving the link to be used; this must
    be one of \code{link = probit} (default option) or \code{link = logit}.}
  \item{knot}{number of knots; by default, \code{knot = 20}.}
  \item{dif}{order of penalty; by default, \code{dif = 2}.}
  \item{raw}{logical; if \code{TRUE}, \code{Y} will be converted to
    unit Fréchet scale. If \code{FALSE}, \code{Y} will be understood as
    already in unit Fréchet scale.}  
}
\value{
  \item{I}{binary vector indicating exceedances.}
  \item{tex}{time of exceedances.}
  \item{traj}{matrix containing posterior simulated trajectories of the
    target of interest.}
  \item{trajhat}{mean trajectory of target of interest.}
  \item{beta}{matrix with posterior samples for the coefficients of B-spline basis.}
  \item{tau2}{matrix with posterior samples for the variance
    parameters.}
  \item{afun}{\code{approxfun} version of fit.}
}
\details{
  Bayesian P-Spline modeling of time-varying extremal dependence was
  introduced by Lee et al (2022). For \code{chit} Gibbs sampler is
  used for \code{link = probit} whereas Metropolis--Hastings with IWLS
  proposal is used for \code{link = logit}. For \code{chibart}, a
  Metropolis--Hastings with IWLS proposal is used for both links. 
}
\references{
  Lee, J., de Carvalho, M., Rua, A., and Avila, J. (2022)
  Bayesian smoothing for time-varying extremal dependence. 
  \emph{Submitted.}
}
\examples{

packages <- c("copula", "evd", "ggplot2")
sapply(packages, require, character.only = TRUE)

n <- 40000 
theta <- function(t)
  (1 - (1 - t)^2) * (cos(4 * pi * (1 - t)^2) + 1) / 4 + 0.2 
t <- seq(0, 1, length.out = n)
XY0 <- data.frame(t, t(sapply(theta(t), function(each) rCopula(1, gumbelCopula(1 / each)))))

fit0 <- chit(XY0, tau = 0.975, link = "probit", prior = list(a = 1, b = 1))
plot(fit0) + 
    geom_line(data = data.frame(x = t, y = 2 - 2^(theta(t))), 
            aes(x = x, y = y), col = "red")


}
