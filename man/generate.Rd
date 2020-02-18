\name{generate}
\alias{generate}

\title{Data simulation}

\description{Generates samples following the Principal Fitted Components model with \eqn{r = 2}, \eqn{d = 1} and covariance matrix {{\eqn{\sigma}}\out{<sup>2</sup>} times the \eqn{p x p} identity matrix (for more details see Bergesio et al. (2020)). The simulated model has covariate
 \eqn{Fy = (y,}\out{y<sup>2</sup>)<sup>t</sup>}}, where \eqn{y} has uniform distribution in the \eqn{(0,4)} interval.}
\usage{
generate(p, n, mutrue, gammatrue, betatrue, sigmatrue)
}

\arguments{
  \item{p}{dimension of the response vector in the inverse model}
  \item{n}{size of the sample required}
  \item{mutrue}{vector of dimension \eqn{p} with the true value of parameter \eqn{\mu} }

  \item{gammatrue}{matrix \eqn{p x d} with the true value of parameter \eqn{\Gamma}}

  \item{betatrue}{matrix \eqn{d x r}, with \eqn{d=1},  containing the true value of the parameter \eqn{\beta}}

\item{sigmatrue}{constant standard deviation of the error term}
}
\value{List with the following components
\item{X}{\eqn{n x p} matrix, each row is a response vector}
\item{Fy}{\eqn{n x r} matrix, each row is the corresponding covariate vector}
}

\references{Bergesio, A., Szretter Noste, M. E. and Yohai, V. J. (2020). A robust proposal of estimation for the
sufficient dimension reduction problem}

\examples{require(MASS)
p=10
generate(p,n=1,mutrue=rep(0,p),gammatrue=as.matrix(c(1,rep(0,p-1))),
betatrue=t(as.matrix(c(1,1))),sigmatrue=1)
}
