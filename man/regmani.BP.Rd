\encoding{UTF-8}
\name{regmani.BP}
\alias{regmani.BP}
%% \alias{regmani.BP.default}
\title{Learning about Regression Manifold from Data}
\description{
  This function computes random Bernstein polynomial estimators of the
  regression manifold of a multivariate extreme value distribution.
}
\usage{
regmani.BP(Y, tau = 0.95, c = 0.0001, m = 'maximize', T = 1000, burn = 100, 
           grid = seq(0.01, 0.99, length = 2^8), raw = TRUE)
}
\arguments{
  \item{Y}{data frame with two columns from which the estimate is to
    be computed.}
  \item{tau}{value used to threshold the data; by default it is set as
    the 0.95 quantile of the pseudo-radius.}
  \item{c}{hyperparameter of Dirichlet prior for weights; see details.}
  \item{m}{string \code{'maximize'} or number of basis functions; see details.}
  \item{T}{number of MCMC iterations; by default \code{T = 1000}.}
  \item{burn}{number of burn in iterations; by default \code{burn = 100}.}
  \item{grid}{grid with coordinates of the points where the angular
    density is estimated; by default \code{grid = seq(0.01, 0.99, length
      = 2^8)}.}
  \item{raw}{logical; if \code{TRUE}, \code{Y} will be converted to
    unit Fréchet scale. If \code{FALSE}, \code{Y} will be understood as
    already in unit Fréchet scale.}
}
\value{
  \item{PI}{matrix containing posterior simulated mixing weights.}
  \item{PIhat}{vector of posterior mean mixing weights.}
  \item{indices}{matrix with indices of Bernstein polynomials (vertex 
    coefficients appear on the first $p$ rows).}
  \item{J}{row sum of Bernstein polynomial index matrix; see details.}
  \item{m}{number of basis functions; see details.}
  \item{w}{pseudo-angles.}
  \item{Y}{raw data.}
  The \code{plot} method depicts a plug-in Bayesian estimator of the
  the angular manifold. Cross sections of the fitted angular manifold
  are plotted when \code{is.numeric(q)} or \code{is.numeric(pred)};
  for example, the conditional median is depicted using \code{q = 0.5},
  and the angular manifold conditional on a specific value of the
  covariate is plotted using \code{pred = a} with \code{a > 0}.
}
\details{
  Regression manifolds were introduced by de Carvalho et al (2022) along 
  with inference procedures based on Hanson et al (2017). 
  The componentwise adaptive Metropolis–Hastings of Haario et al (2005) is used 
  to fit the model. Specifically, a static MH update is used over the burn-in 
  period, between \code{burn} and \code{2 * burn} the algorithm starts
  building adaptive variance, and after \code{2 * burn} it uses adaptive
  variance and it harvests iterates. This implies that posterior means
  are based on \code{T - 2 burn} iterates, and hence \code{burn} and
  \code{T} need to be set so that \code{T > 2 burn}. A Dirichlet prior,
  Dir(\eqn{c, \dots, c}{c,...,c}), is assigned to the Bernstein polynomial 
  weights, and following Hanson et al (2017), when \code{m =
    'maximize'},  \eqn{J}{J} is chosen so be maximized 
  s.t. \eqn{m \leq k}{m <= k}, with \eqn{k}{k} denoting the number of
  pseudo-angles.   
}
\references{
  de Carvalho, M., Kumukova, A. and dos Reis, G. (2022) Regression-type analysis 
  for multivariate extreme values. Submitted.
  
  Haario, H., Eero S. and Johanna T. (2005) Componentwise adaptation for    
  high dimensional MCMC. \emph{Computational Statistics}, 20, 265--73.
  
  Hanson, T., de Carvalho, M. and Chen, Y. (2017) Bernstein polynomial      
  angular densities of multivariate extreme value distributions.            
  \emph{Statistics and Probability Letters}, 128, 60--66.
}
\examples{
## ================= ## 
##  Simulated data   ##
## ================= ## 
set.seed(1)
n <- 5000
alpha <- 10 # = 1 / lambda
ones <- c(1, 1, 1)
Y <- evd::rbvevd(n, model = "hr", mar1 = ones, mar2 = ones, dep = alpha)
fit <- regmani.BP(Y)

## plot regression manifold
pfit <- plot(fit, main = "Estimate", zlim = c(0, 45))
ptrue <- plot(regmani("hr", dep = alpha), main = "True", zlim = c(0, 45))
print(pfit, split = c(1, 1, 2, 1), more = TRUE)
print(ptrue, split = c(2, 1, 2, 1), more = FALSE)

## plot cross sections of regression manifold
plot(fit, q = c(0.1, 0.5, 0.9))
plot(fit, pred = c(1, 5, 10))

## plot cross sections of regression manifold against true
true <- regmani("hr", dep = alpha, q = c(0.1, 0.5, 0.9))
plot(fit, q = c(0.1, 0.5, 0.9)) + 
  ggplot2::geom_line(data = true$rm, ggplot2::aes(x = x, y = z),
                     linetype = "dashed")

\dontrun{
## plot cross sections along with credible bands
plot(fit, q = c(0.1, 0.5, 0.9), bands = TRUE)
plot(fit, x = c(1, 5, 10), bands = TRUE)}
}
