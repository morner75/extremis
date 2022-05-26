\encoding{UTF-8}
\name{cPTdensity}
\alias{scedasis.PT}
%% \alias{scedasis.PT.default}
\title{Nonparametric Bayesian Estimation of Structure Scedasis Density using Mixtures of Polya Trees}
\description{
  This function computes a Bayesian structure scedasis density estimate generating a posterior sample of a finite
    Mixture of Polya trees.% For the univariate case it computes Einmahl et al (2016) approach using a finite mixture of Polya trees.
}
\usage{
scedasis.PT(XY, tau = 0.95, raw = TRUE, structure = "min", prior, mcmc, state,
 status, data = sys.frame(sys.parent()), na.action = na.fail)
}

\arguments{
 \item{XY}{data frame from which the estimate is to be computed; first
   column corresponds to time, the second and third columns correspond
   to the variables of interest \code{X} and \code{Y} respectively.}
 \item{tau}{value used to threshold the data \code{z = min(X, Y)}; by default
   \code{threshold = quantile(z, tau)}.}
 \item{raw}{logical; if \code{TRUE}, \code{X} and \code{Y} will be
   converted to unit Fréchet scale. If \code{FALSE}, \code{X} and
   \code{Y} will be understood as already in unit Fréchet scale. If a
   single variable is provided
   \code{raw} is automatically set to \code{TRUE}.}
 \item{structure}{the structure scedasis option: "min","max","sum". The default will compute the minimum.}
 \item{prior}{a list giving the prior information. The list includes the
   following parameter: \code{a0} and \code{b0} giving the
   hyperparameters for prior distribution of the precision parameter of
   the Poly tree prior, \code{alpha} giving the value of the precision
   parameter (it must be specified if \code{alpha} is missing, see
   details below), optionally \code{M} giving the finite level to be
   considered (if \code{M} is specified, a partially specified mixture
   of Polya trees model is fitted), \code{tau1} and \code{tau2} giving
   the hyperparameters of the inverted gamma prior distribution for the
   centering Beta scale parameter,\code{m0} and \code{S0} giving the
   hyperparameters of the inverted gamma prior distribution for the
   centering Beta scale parameter, and \code{al}, \code{be} giving the
   value of the shape and scale parameter of the centering
   distribution.}
 \item{mcmc}{a list giving the MCMC parameters. The list must include
   the following integers: \code{nburn} giving the number of burn-in
   scans, \code{nskip} giving the thinning interval, \code{nsave} giving
   the total number of scans to be saved, \code{ndisplay} giving the
   number of saved scans to be displayed on screen (the function reports
   on the screen when every \code{ndisplay} iterations have been carried
   out), \code{tune1}, \code{tune2}, and \code{tune3} giving the
   positive Metropolis tuning parameter for the baseline shape, scale,
   and precision parameter, respectively (the default value is 1.1)}
\item{state}{a list giving the current value of the parameters. This
  list is used if the current analysis is the continuation of a previous
  analysis.}
\item{status}{a logical variable indicating whether this run is new
  (\code{TRUE}) or the continuation of a previous analysis
  (\code{FALSE}). In the latter case the current value of the parameters
  must be specified in the object \code{state}.}
\item{data}{data frame.}
\item{na.action}{a function that indicates what should happen when the data
                     contain \code{NA}s. The default action (\code{na.fail}) causes
                     \code{PTdensity} to print an error message and terminate if there are any
                     incomplete observations.}
}

\value{
  \item{c}{Structure scedasis density estimator.}
  \item{k}{number of exceedances above the threshold.}
  \item{w}{standardized indices of exceedances.}
  \item{Y}{raw data.}
  \item{al}{ giving the value of the baseline shape parameter.}
  \item{be}{ giving the value of the baseline scale parameter.}
  \item{alpha}{ giving the value of the precision parameter.}
  The \code{plot} method depicts the estimated structure scedasis density.
}
\details{This function learns about the structure scedasis function using a
  Mixture of Polya Trees (MPT) prior as proposed in Palacios and de Carvalho (2020). In the
  particular case where XY contains no third column, the function learns
  about the scedasis function of Einmahl et al (2016) using an MPT
  prior. The details are as follows. Let \deqn{Z_i = min(X_i, Y_i)}.
  The model is given by:
  \deqn{Z_1, \ldots , Z_n | G \sim G}{Z1,...,Zn | G ~ G}
  \deqn{G | \alpha,al,be \sim PT(\Pi^{al,be},\textit{A})}{G | alpha,a,b ~ PT(Pi^{al,be},\textit{A})}
  where, the the PT is centered around a \eqn{Beta(al,be)}{Beta(al,be)} distribution, by taking each \eqn{m} level of the partition \eqn{\Pi^{al, be}}{\Pi^{al, be}} to coincide
  with the \eqn{k/2^m, k=0,\ldots,2^m} quantile of the \eqn{Beta(al,be)}{Beta(al,be)} distribution.
  The family \eqn{\textit{A}=\{\alpha_e: e \in E^{*}\}}{\textit{A}=\{alphae: e \in E^{*}\}},
  where \eqn{E^{*}=\bigcup_{m=0}^{m} E^m}{E^{*}=\bigcup_{m=0}^{M} E^m}
  and \eqn{E^m} is the \eqn{m}-fold product of \eqn{E=\{0,1\}},
  was specified as \eqn{\alpha_{e_1 \ldots e_m}=\alpha m^2}{alpha{e1 \ldots em}=\alpha m^2}.

  In the univariate case,the following proper priors can be assigned:
    \deqn{al |  m_0, S_0 \sim LNormal(m_0,S_0)}{al | m0, S0 ~ LNormal(m0,S0)}
  \deqn{be | \tau_1, \tau_2 \sim LNormal(\tau_1,\tau_2)}{be | tau1, tau2 ~ LNormal(tau1,tau2)}

   To complete the model specification, independent hyperpriors are assumed,
  \deqn{\alpha | a_0, b_0 \sim \Gamma(a_0,b_0)}{alpha | a0, b0 ~ \Gamma(a0,b0)}

  The precision parameter, \eqn{\alpha}{alpha}, of the \code{PT} prior
  can be considered as random, having a \code{gamma} distribution, \eqn{\Gamma(a_0,b_0)}{\Gamma(a0,b0)},
  or fixed at some particular value. To let \eqn{\alpha}{alpha} to be fixed at a particular
  value, set \eqn{a_0}{a0} to NULL in the prior specification.

  In the computational implementation of the model, Metropolis--Hastings steps are used to
  sample the posterior distribution of the baseline and precision parameters.
}

\references{
  Palacios, V., de Carvalho, M. (2020) Bayesian semiparametric modeling
  of jointly heteroscedastic extremes. \emph{Preprint.}
  %% Einmahl, J. H., Haan, L., and Zhou, C. (2016) Statistics of
  %% heteroscedastic extremes. \emph{Journal of the Royal Statistical
  %% Society: Ser. B}, 78(1), 31--51.

  %% Hanson, T. (2006) Inference for Mixtures of Finite Polya Trees.
  %% \emph{ Journal of the American Statistical Association}, 101: 1548-1565.
 % Lavine, M. (1992) Some aspects of Polya tree distributions for statistical
  %modelling. The Annals of Statistics, 20: 1222-11235.
  %Lavine, M. (1994) More aspects of Polya tree distributions for statistical
  %modelling. The Annals of Statistics, 22: 1161-1176.
  }

\examples{
\dontrun{
## Example from Palacios and de Carvalho (2020, submitted)
library(evd)
## Initial state
state <- NULL
## MCMC parameters
nburn <- 2000
nsave <- 1000
nskip <- 0
ndisplay <- 500
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay,
 tune1=1.1,tune2=1.1,tune3=1.1)
## Prior information
prior<-list(a0=1,b0=1,M=8,m0=.01,S0=.01,tau1=.01,tau2=.01);

T <- 5000
time <- seq(1/T, 1, by = 1/T)
set.seed(12333)
aux <- matrix(0, T, 2)
for (i in 1:T) {
    aux[i,]<-rbvevd(1, dep =.5, asy=c(sin(pi*time[i]),time[i]),
    model="alog",mar1 = c(1, 1, 1), mar2 = c(1, 1, 1))
}
XY <- cbind(time, aux)
fit <- scedasis.PT(XY, prior = prior, mcmc = mcmc, state = state, status =
TRUE)
plot(fit)
}}
