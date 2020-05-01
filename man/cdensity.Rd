\encoding{UTF-8}
\name{cdensity}
\alias{cdensity}
\alias{cdensity.default}
\title{Kernel Smoothed Scedasis Density}
\description{
  This function computes a kernel scedasis density estimate.
}
\usage{
cdensity(Y, threshold = quantile(Y[, 2], 0.95), ...)
}

\arguments{
  \item{Y}{data frame from which the estimate is to be computed; first
    column corresponds to time and the second to the variable of interest.}  
  \item{threshold}{value used to threshold the data \code{y}; by default
    \code{threshold = quantile(y, 0.95)}.}
  \item{...}{further arguments for \code{density} methods.}
}

\value{
  \item{c}{scedasis density estimator.}
  \item{k}{number of exceedances above the threshold.}
  \item{w}{standardized indices of exceedances.}
  \item{Y}{raw data.}
  The \code{plot} method depicts the smooth scedasis density. 
}
\details{
  Kernel smoothing for the scedasis density was introduced by
  Einmahl et al (2016). %% The option \code{unit = "true"} was suggested by de
  %% Carvalho (2016), as a way to prevent putting mass outside the unit
  %% interval.
}

\references{
  %% de Carvalho, M., (2016)
  %% Statistics of extremes: Challenges and opportunities. In: \emph{Handbook of
  %% EVT and its Applications to Finance and Insurance} Eds
  %% F. Longin. Hoboken: Wiley.

  Einmahl, J. H., Haan, L., and Zhou, C. (2016) Statistics of
  heteroscedastic extremes. \emph{Journal of the Royal Statistical
  Society: Ser. B}, 78(1), 31--51.
}
\examples{
data(lse)
attach(lse)
Y <- data.frame(DATE[-1], -diff(log(ROYAL.DUTCH.SHELL.B)))
T <- dim(Y)[1]
k <- floor((0.4258597) * T / (log(T)))
fit <- cdensity(Y, kernel = "biweight", bw = 0.1 / sqrt(7), 
                threshold = sort(Y[, 2])[T - k])
plot(fit)
plot(fit, original = FALSE)
}
\author{Miguel de Carvalho}