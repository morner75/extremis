\encoding{UTF-8}
\name{kgvar}
\alias{kgvar}
\alias{kgvar.default}
\title{K-Geometric Means Algorithm for Value-at-Risk}
\description{
  This function performs k-geometric means for time-varying
  value-at-risk.
}
\usage{
  kgvar(y, centers, iter.max = 10, conf.level = 0.95)
}
\arguments{
  \item{y}{data frame from which the estimate is to be computed; first
    column corresponds to time and the second to the remainder of
    interest.}
  \item{centers}{the number of clusters or a set of initial
    (distinct) cluster centres. If a number, a random set of (distinct)
    rows in \code{y} is chosen as the initial centers.}
   \item{iter.max}{the maximum number of iterations allowed. The
   default is 10.}
   \item{conf.level}{the confidence level. The default is 0.95.}
}
\value{kgvar returns an object of class \code{"kgvar"} which has a
  fitted method.  It is a list with at least the following components:
  \item{var.new}{cluster center value-at-risk function.}
  \item{clusters}{cluster allocation.}
  \item{Y}{raw data.}
  \item{n.clust}{number of clusters.}
  \item{scale.param}{the scale parameters in the Pareto-like tail
    specification.}
    \item{conf.level}{the confidence level.}
    \item{hill}{hill estimator of extreme value index.}
  The \code{plot} method depicts the k-geometric means algorithm for
  time-varying value-at-risk. If \code{c.c} is \code{TRUE}, the method displays the
  cluster means.
}
\details{The intermediate sequence \eqn{\kappa_T} is chosen
  proportional to \eqn{T/\log T}.
}
\references{Rubio, R., de Carvalho, M. and Huser, R. (2018)
  Similarity-Based Clustering of Extreme Losses from the London Stock
  Exchange. Submitted.
}
\examples{
\dontrun{
## Example (Overlapping version of Fig. 8 in Supplementary Materials)
data(lse)
attach(lse)
y <- -apply(log(lse[, -1]), 2, diff) 
fit <- kgvar(y, centers = 3)
plot(fit, c.c = TRUE, ylim = c(0, 0.1))
}
}
\author{Miguel de Carvalho, Rodrigo Rubio.}