\name{cross_val_pfc}
\alias{cross_val_pfc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selection of the dimension of the sufficient reduction subspace
%%  ~~function to do ... ~~
}
\description{Computes the dimension of the sufficient reduction subspace by cross-validation
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cross_val_pfc(X, Fy, dmax, aux, grafico=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{X}{vector of response variables in the inverse model
%%     ~~Describe \code{x} here~~
}
\item{Fy}{vector of covariates in the inverse problem, vector containig functions of the response variable in the original problem
%%
}
  \item{dmax}{maximum value of the dimension to be evaluated
%%     ~~Describe \code{dmax} here~~
}
  \item{aux}{list containing the constants for the tau-scale
%%     ~~Describe \code{aux} here~~
}

  \item{grafico}{if TRUE, shows graphs of the Cross-Validation function versus dimension
%%     ~~Describe \code{grafico} here~~
}
}
\details{
The dimension of the reduction subspace, using both maximum likelihood
and robust tau estimators to select is computed using five fold cross validation for the PFC model. See Bergesio et al. (2020).

Computes \code{sd} of both objective functions and gives back the smaller value that satisfies the
one standard deviation rule from the minimum. See Hastie et al. (2009)
%%  ~~ If necessary, more details than the description above ~~
}
\value{List with the following components

\item{d.crossval.ml }{dimension selected by cross-validation using ML estimator}
\item{obj.ml}{objective function values for the ML estimator, for values of d between 0 and \code{dmax}}
\item{obj.ml.sd}{standard deviation of the objective function for the ML estimator, for values of d between 0 and \code{dmax}}
\item{d.crossval.rob}{dimension selected by cross-validation using tau estimator}
\item{obj.rob}{objective funcion values for the robust estimator}
\item{obj.rob.sd}{standard deviation of the objective funcion for the robust estimator}
...
}
\references{
Hastie, T., Tibshirani, R., Friedman, J. H. (2009). The elements of statistical learning: data mining, inference, and prediction. 2nd ed. New York: Springer.

Bergesio, A., Szretter Noste, M. E. and Yohai, V. J. (2020). A robust proposal of estimation for the
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
cross_val_pfc(X, Fy, dmax=4, aux, grafico=TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
