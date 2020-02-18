\name{rhoc}
\alias{rhoc}
\title{Quasi-optimal rho function
}
\description{
rhoc computes the \eqn{\rho} quasi-optimal function defined in the references
}
\usage{
rhoc(x, c)
}
\arguments{
  \item{x}{value where the rho function is evaluated
}
  \item{c}{tuning parameter, it defines where the increasing function starts to be constant
}
}
\references{
Muler,  N. and Yohai,  V. J. (2002). Robust estimates for arch processes.
{\emph{Journal of Time Series Analysis}}, 23(3):341-375.

Salibian-Barrera, M., Van Aelst, S. and Yohai, V. J. (2016). Robust tests for linear regression models based on \eqn{\tau}-estimates. {\emph{Computational Statistics & Data Analysis}} 93:436-455.

}
\examples{
plot(seq(0,1.3,.01),rhoc(seq(0,1.3,.01),1),type='l',ylab="",xlab="")
}
