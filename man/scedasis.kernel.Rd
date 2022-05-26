\encoding{UTF-8}
\name{cdensity.kernel}
\alias{scedasis.kernel}
%% \alias{scedasis.kernel.default}
\title{Kernel Smoothed Joint Scedasis Density}
\description{
  This function computes a kernel joint scedasis density estimate.
}
\usage{
scedasis.kernel(XY, tau = 0.95, raw = TRUE, structure = "min", ...)
}

\arguments{
 \item{XY}{data frame from which the estimate is to be computed; first
   column corresponds to time, the second and third columns correspond
   to the variables of interest \code{X} and \code{Y} respectively.}
 \item{tau}{value used to threshold the data \code{Z = min(X, Y)}; by default
   \code{threshold = quantile(Z, tau)}.}
 \item{raw}{logical; if \code{TRUE}, \code{X} and \code{Y} will be
   converted to unit Fréchet scale. If \code{FALSE}, \code{X} and
   \code{Y} will be understood as already in unit Fréchet scale. If a
   single variable is provided
   \code{raw} is automatically set to \code{TRUE}.}  
 \item{structure}{the structure scedasis option: "min","max","sum".
   The default will compute the minimum.}
  \item{...}{further arguments for \code{density} methods.}  
}

\value{
  \item{c}{scedasis density estimator.}
  \item{k}{number of exceedances above the threshold.}
  \item{w}{standardized indices of exceedances.}
  \item{XY}{raw data.}
  The \code{plot} method depicts the smooth scedasis density. 
}

\details{
  This function learns about the joint scedasis function using kernel
  methods as discussed in Palacios and de Carvalho (2022). In the
  particular case where \code{XY} contains no third column, the function learns
  about the scedasis function of Einmahl et al (2016). 
  %% The option \code{unit = "true"} was suggested by de
  %% Carvalho (2016), as a way to prevent putting mass outside the unit
  %% interval.
}
\references{
  de Carvalho, M., (2016)
  Statistics of extremes: Challenges and opportunities. In: \emph{Handbook of
    EVT and its Applications to Finance and Insurance} Eds F. Longin.
  Hoboken: Wiley.
  
  Einmahl, J. H., Haan, L., and Zhou, C. (2016) Statistics of
  heteroscedastic extremes. \emph{Journal of the Royal Statistical
    Society: Ser. B}, 78(1), 31--51.
  
  Palacios, V. and de Carvalho, M. (2022) Bayesian semiparametric
  modeling of jointly heteroscedastic extremes. \emph{Preprint.}
}
\examples{
data(lse)
attach(lse)
XY <- data.frame(DATE[-1], -diff(log(ROYAL.DUTCH.SHELL.B)))
T <- dim(XY)[1]
k <- floor((0.4258597) * T / (log(T)))
fit <- scedasis.kernel(XY, kernel = "biweight", bw = 0.1 / sqrt(7))
plot(fit)
plot(fit, original = FALSE)

## Example from Palacios and de Carvalho (2022, submitted)
library(evd)
T <- 5000
time <- seq(1/T, 1, by = 1/T)
set.seed(1263)
aux <- matrix(0,T,2)
for(i in 1:T) {
    aux[i,] <- rbvevd(1,dep=sin(time[i]*pi),model="log",
    mar1 = c(1, 1, 1), mar2 = c(1,1,1))}
XY <- cbind(time,aux)
fit <- scedasis.kernel(XY, kernel = "biweight", bw = 0.1)
plot(fit)
}
