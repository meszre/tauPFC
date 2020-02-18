\name{psic}
\alias{psic}
\title{First derivative of the quasi-optimal rho function
}
\description{
psic computes the first derivative of the \eqn{\rho} quasi-optimal function defined in the references
}
\usage{
psic(x, c)
}
\arguments{
  \item{x}{value where the psic function is evaluated
}
  \item{c}{tuning parameter, it defines where the \eqn{\rho} quasi-optimal function starts to be constant
}
}
\references{
Muler,  N. and Yohai,  V. J. (2002). Robust estimates for arch processes.
{\emph{Journal of Time Series Analysis}}, 23(3):341-375.

Salibian-Barrera, M., Van Aelst, S. and Yohai, V. J. (2016). Robust tests for linear regression models based on \eqn{\tau}-estimates. {\emph{Computational Statistics & Data Analysis}} 93:436-455.
}
\examples{
plot(seq(0,1.3,.01),psic(seq(0,1.3,.01),1),type='l',ylab="",xlab="")
}
