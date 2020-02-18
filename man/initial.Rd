\name{initial}
\alias{initial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initial coefficients and covariance estimators for PFC model
%%  ~~function to do ... ~~
}
\description{
Compute initial estimators for the coefficient matrix and the covariance matrix. Both are strongly robust to the presence of outliers in the sample but possibly inefficient. They are a good starting point for tau-estimators for PFC model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
initial(X, Fy, aux, efficiency = 0.85)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{vector of response variables in the inverse model, {\eqn{n x p}} matrix, each row is a response vector
%%     ~~Describe \code{x} here~~
}
\item{Fy}{vector of covariates in the inverse problem, vector containing functions of the response variable in the original multiple regression problem. Is a {\eqn{n x p}} matrix, each row is the corresponding response vector
%%
}
  \item{aux}{list containing the constants for the tau-scale (with components c1, k1, c2, k2, as described in  %\code{\link[kappa_and_c]}{kappa_and_c}
  \code{\link[tauPFC]}{kappa_and_c})
  %\href{kappa_and_c}{kappa_and_c})
%%     ~~Describe \code{aux} here~~
}
  \item{efficiency}{required efficiency for initial robust estimator
%%     ~~Describe \code{efficiency} here~~
}
}
\details{This function computes the coefficient matrix concatenating the coefficients resulting from the fitting of univariate regressions using a robust estimator with high breakdown point and high efficiency (require "robust"). Univariate estimators are computed using  %\href{robust::lmRob.Rd}{lmRob} \code{\link[robust]}{lmRob}
"lmRob" function. It computes a robust estimator of errors covariance matrix (require "rrcov").
%  ~~ If necessary, more details than the description above ~~
}
\value{
The estimation is given as a list with components
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{beta0 }{coefficient matrix }
\item{delta0 }{S covariance matrix of residuals}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{For the final estimation use \code{\link{tauestimate}} with this initial estimator as an input
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
p=10
mutrue=rep(0,p)
gamatrue=as.matrix(c(1,rep(0,p-1)))
betatrue=t(as.matrix(c(1,1)))

data_sim=generate(p,n=200,mutrue,gamatrue,betatrue,sigmatrue=1)

aux=kappa_and_c(p,efficiency=0.85)
initial(data_sim$X,data_sim$Fy,aux,efficiency=0.85)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
