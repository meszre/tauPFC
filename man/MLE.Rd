\name{MLE}
\alias{MLE}

\title{Maximum likelihood estimators for the Principal Fitted Components (PFC) Model}
\description{Computes the ML estimators for the PFC model
}
\usage{
MLE(X,Fy,d)
}
\arguments{
\item{X}{vector of response variables in the inverse model, {\eqn{n x p}} matrix, each row is a response vector}
\item{Fy}{vector of covariates in the inverse problem, vector containing functions of the response variable in the original multiple regression problem. Is a {\eqn{n x r}} matrix, each row is the corresponding response vector
}

\item{d}{number indicating the reduction subspace dimension
}
}
\details{We consider the Principal Fitted  Components (PFC) model given by
\eqn{X = \mu + \Gamma\beta f(y) +\Delta}\out{<sup>1/2</sup>}\eqn{\epsilon},
where the variables are
\itemize{
\item{} \eqn{y} is an observed response variable of the original model. \eqn{f} is a known vector valued function, that takes values in R\out{<sup>r</sup>}
\item{} \eqn{X} is the correspondent \emph{p x 1} observed covariates vector
\item \eqn{\epsilon} is unobserved  \emph{p} dimensional vector, \eqn{\Delta}\out{<sup>1/2</sup>}\eqn{\epsilon} is the error vector}
and the unknown parameters (to be estimated) are
\itemize{
\item \eqn{\mu} a \emph{p x 1} vector of intercepts
\item{} \eqn{\Gamma} is a full-rank \emph{p x d} matrix whose columns span the dimension reduction subspace
\item{} \eqn{\beta} is a full-rank \emph{d x r} matrix
\item{} cov(\eqn{\epsilon}) = \eqn{\Delta}, is a \emph{p x p} positive definite  matrix }
Both coefficient matrices \eqn{\Gamma} and \eqn{\beta} are not unique, but their product \emph{p x r} matrix is unique, with rank \code{d} \eqn{\le min(p,r)}. The notation refers to Cook and Forzani (2008).}

\value{
List with the following components
\item{mu }{MLE estimation of the term \eqn{\mu} in the PFC model}
\item{beta }{MLE estimation of the parameter \eqn{\beta} in the PFC model}
\item{gamma }{MLE estimation of the parameter \eqn{\Gamma} in the PFC model}
\item{delta}{MLE estimation of the covariance matrix \eqn{\Delta} in the PFC model}
}
\references{
Cook, R. D. and Forzani, L. (2008). Principal Fitted components for dimension reduction in
regression. Statistical Science, 23(4):485-501.
}
\examples{
p=10
n=200
mutrue=rep(0,p)
gamatrue=as.matrix(c(1,rep(0,p-1)))
betatrue=t(as.matrix(c(1,1)))

data_sim=generate(p,n,mutrue,gamatrue,betatrue,sigmatrue=1)
Fy=data_sim$Fy
X=data_sim$X
MLE(X,Fy,d=1)

}
