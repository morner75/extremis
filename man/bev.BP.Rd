\encoding{UTF-8}
\name{bev.BP}
\alias{bev.BP}
%% \alias{bev.BP.default}
\title{Random Bernstein Polynomial Inference for Multivariate Extremes}
\description{
  This function computes random Bernstein polynomials estimators of the
  angular density and angular distribution function of a multivariate
  extreme value distribution.
}
\usage{
bev.BP(Y, tau = 0.95, c = 0.0001, m = 'maximize', T = 1000, burn = 100, 
       grid = seq(0.01, 0.99, length = 2^8), target = "pdf", raw = TRUE)
}
\arguments{
  \item{Y}{data frame with two columns from which the estimate is to
    be computed.}
  \item{tau}{value used to threshold the data; by default it is set as
    the 0.95 quantile of the pseudo-radius.}
  \item{c}{hyperparameter of Dirichlet prior for weights; see details}
  \item{nu}{concentration parameter of beta distribution which
    controls the amount of smoothing.}
  \item{T}{number of MCMC iterations; by default \code{T = 1000}.}
  \item{burn}{number of burn in iterations; by default \code{burn = 100}.}
  \item{grid}{grid with coordinates of the points where the angular
    density is estimated; by default \code{grid = seq(0.01, 0.99, length
      = 2^8)}.}
  \item{target}{value set to \code{"pdf"} or \code{"cdf"} depending on
    whether the goal is to estimate an angular density or an angular
    distribution function, respectively; by default: \code{target =
      "pdf"}.}
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
  \item{traj}{matrix containing posterior simulated trajectories of the
    angular density (\code{target = "pdf"}), distribution function
    (\code{target = "cdf"}) computed over \code{grid}.}
  \item{trajhat}{mean trajectory of target of interest evaluated at
    \code{grid}.}
  \item{w}{pseudo-angles.}
  \item{target}{string with target of interest.}
  The \code{plot} method depicts the estimated object [angular density
  (\code{target = 'pdf'}) or angular distribution function (\code{target =
    'cdf'})] along with 95\% probability bands, unless \code{bands =
    FALSE}. 
}
\details{
  Random bernstein angular densities were introduced by Hanson et al (2017). 
  The componentwise adaptive Metropolis–Hastings of Haario et al (2005) is used 
  to fit the model. Specifically, a static MH update is used over the burn-in 
  period, between \code{burn} and \code{2 * burn} the algorithm starts
  building adaptive variance, and after \code{2 * burn} it uses adaptive
  variance and it harvests iterates. This implies that posterior means
  are based on \code{T - 2 burn} iterates, and hence \code{burn} and
  \code{T} need to be set so that \code{T > 2 burn}. A Dirichlet prior,
  Dir(\eqn{c, \dots, c}{c,...,c}), is assigned to the Bernstein polynomial 
  weights, and following Hanson et al (2017), when \code{m =
    'maximize'}, \eqn{J}{J} is chosen so be maximized 
  s.t. \eqn{m \leq k}{m <= k}, with \eqn{k}{k} denoting the number of
  pseudo-angles.   
}
\references{
  Haario, H., Eero, S. and Johanna, T. (2005) Componentwise adaptation for    
  high dimensional MCMC. \emph{Computational Statistics}, 20, 265--73.
                                                                               
  Hanson, T., de Carvalho, M. and Chen, Y. (2017) Bernstein polynomial      
  angular densities of multivariate extreme value distributions.            
  \emph{Statistics and Probability Letters}, 128, 60--66.
}
\examples{
## ================= ## 
##  Beatenberg data  ##
## ================= ## 
data(beatenberg)
fit <- bev.BP(beatenberg, tau = 0.98, raw = FALSE)
plot(fit, target = "pdf")
plot(fit, target = "cdf")
\dontrun{
## Run the lines below to compare with the kernel angular density
## estimator of de Carvalho et al (2013, Fig. 7)         
fit.kernel <- bev.kernel(beatenberg, tau = 0.98, nu = 163, raw = FALSE)
chart <- plot(fit) +
  ggplot2::geom_line(ggplot2::aes(x = fit.kernel$grid,
                                  y = fit.kernel$estimate),
                                  lty = 2, col = "red") 
plot(chart)

## Run the lines below to compare with the kernel angular measure 
## estimator of de Carvalho et al (2013, Fig. 7)
fit <- bev.BP(beatenberg, tau = 0.98, raw = FALSE, target = "cdf")
fit.kernel <- bev.kernel(beatenberg, tau = 0.98, nu = 163,
                         raw = FALSE, target = "cdf")
chart <- plot(fit) + 
  ggplot2::geom_line(ggplot2::aes(x = fit.kernel$grid,
                                  y = fit.kernel$estimate),
                                  lty = 2, col = "red") 
plot(chart)
}
\dontrun{
## ================= ## 
##  Simulated data   ##
## ================= ## 
n <- 5000
alpha <- 10 # = 1 / lambda
ones <- c(1, 1, 1)
Y <- evd::rbvevd(n, model = "hr", mar1 = ones, mar2 = ones, dep = alpha)
fit <- bev.BP(Y, tau = 0.98)
chart <- plot(fit)
## now let's compare with the true angular density 
grid <- seq(0.01, 0.99, by = 0.01)
true <- function(x)
  evd::hbvevd(x, model = "hr", dep = alpha) / 2
plot(fit) +
  ggplot2::geom_line(ggplot2::aes(x = grid, y = true(grid)), lty = 2, col = "red") 
}
}
