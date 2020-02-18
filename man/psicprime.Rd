\name{psicprime}
\alias{psicprime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Second derivative of the quasi-optimal rho function
%%  ~~function to do ... ~~
}
\description{psicprime computes the second derivative of the \eqn{\rho} quasi-optimal function defined in the references
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
psicprime(x, c)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{value where the second derivative function is evaluated
%%     ~~Describe \code{d} here~~
}
  \item{c}{tuning parameter, it defines where the \eqn{\rho} quasi-optimal function starts to be constant
%%     ~~Describe \code{c} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Muler,  N. and Yohai,  V. J. (2002). Robust estimates for arch processes.
{\emph{Journal of Time Series Analysis}}, 23(3):341-375.

Salibian-Barrera, M., Van Aelst, S. and Yohai, V. J. (2016). Robust tests for linear regression models based on \eqn{\tau}-estimates. {\emph{Computational Statistics & Data Analysis}} 93:436-455.

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
\examples{
plot(seq(0,1.3,.01),psicprime(seq(0,1.3,.01),1),type='l',ylab="",xlab="")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
