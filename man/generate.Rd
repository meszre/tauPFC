\name{generate}
\alias{generate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data simulation
%%  ~~function to do ... ~~
}
\description{Generates samples following the Principal Fitted Components model with \eqn{r = 2}, \eqn{d = 1} and covariance matrix {{\eqn{\sigma}}\out{<sup>2</sup>} times the \eqn{p x p} identity matrix (for more details see Bergesio et al. (2020)). The simulated model has covariate
 \eqn{Fy = (y,}\out{y<sup>2</sup>)<sup>t</sup>}}, where \eqn{y} has uniform distribution in the \eqn{(0,4)} interval.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generate(p, n, mutrue, gammatrue, betatrue, sigmatrue)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{dimension of the response vector in the inverse model
%%     ~~Describe \code{p} here~~
}
  \item{n}{size of the sample required
%%     ~~Describe \code{n} here~~
}
  \item{mutrue}{vector of dimension \eqn{p} with the true value of parameter \eqn{\mu} }

  \item{gammatrue}{matrix \eqn{p x d} with the true value of parameter \eqn{\Gamma}
%%     ~~Describe \code{gammatrue} here~~
}
  \item{betatrue}{matrix \eqn{d x r}, with \eqn{d=1},  containing the true value of the parameter \eqn{\beta}
%%     ~~Describe \code{betatrue} here~~
}
\item{sigmatrue}{constant standard deviation of the error term
%%     ~~Describe \code{betatrue} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{List with the following components
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{X}{\eqn{n x p} matrix, each row is a response vector}
\item{Fy}{\eqn{n x r} matrix, each row is the corresponding covariate vector}
%% ...
}
\references{Bergesio, A., Szretter Noste, M. E. and Yohai, V. J. (2020). A robust proposal of estimation for the
sufficient reduction problem
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{require(MASS)
p=10
generate(p,n=1,mutrue=rep(0,p),gammatrue=as.matrix(c(1,rep(0,p-1))),betatrue=t(as.matrix(c(1,1))),sigmatrue=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
