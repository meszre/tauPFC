\name{cross_val_pfc}
\alias{cross_val_pfc}
\title{
Selection of the dimension of the sufficient reduction subspace}
\description{Computes the dimension of the sufficient reduction subspace by cross-validation}
\usage{
cross_val_pfc(X, Fy, dmax, aux, grafico = TRUE)
}

\arguments{
 \item{X}{vector of response variables in the inverse model}
\item{Fy}{vector of covariates in the inverse problem, vector containing functions of the response variable in the original problem}
  \item{dmax}{maximum value of the dimension to be evaluated}
  \item{aux}{list containing the constants for the tau-scale}

  \item{grafico}{a logical, if TRUE, shows graphs of the Cross-Validation function versus dimension}
}

\details{
The dimension of the reduction subspace, using both maximum likelihood
and robust tau estimators to select is computed using five fold cross validation for the PFC model. See Bergesio et al. (2020).

Computes standard deviations of both objective functions and gives back the smaller value that satisfies the
one standard deviation rule from the minimum. See Hastie et al. (2009)}
\value{List with the following components

\item{d.crossval.ml }{dimension selected by cross-validation using ML estimator}
\item{obj.ml}{objective function values for the ML estimator, for values of d between 0 and \code{dmax}}
\item{obj.ml.sd}{standard deviation of the objective function for the ML estimator, for values of d between 0 and \code{dmax}}
\item{d.crossval.rob}{dimension selected by cross-validation using tau estimator}
\item{obj.rob}{objective function values for the robust estimator}
\item{obj.rob.sd}{standard deviation of the objective function for the robust estimator}
...
}
\references{
Hastie, T., Tibshirani, R., Friedman, J. H. (2009). The elements of statistical learning: data mining, inference, and prediction. 2nd. ed. New York: Springer.

Bergesio, A., Szretter Noste, M. E. and Yohai, V. J. (2020). A robust proposal of estimation for the
sufficient dimension reduction problem}

\examples{
p=10
r=2
n=200
mutrue=rep(0,p)
gamatrue=as.matrix(c(1,rep(0,p-1)))
betatrue=t(as.matrix(c(1,1)))

data_sim=generate(p,n,mutrue,gamatrue,betatrue,sigmatrue=1)
Fy=data_sim$Fy
X=data_sim$X

aux=kappa_and_c(p,efficiency=.9)
cross_val_pfc(X, Fy, dmax=3, aux, grafico=TRUE)
}
