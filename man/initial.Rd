\name{initial}
\alias{initial}
\title{Initial coefficients and covariance estimators for PFC model
}
\description{
Compute initial estimators for the coefficient matrix and the covariance matrix. Both are strongly robust to the presence of outliers in the sample but possibly inefficient. They are a good starting point for tau-estimators for PFC model.
}
\usage{
initial(X, Fy, aux, efficiency = 0.85)
}
\arguments{
  \item{X}{vector of response variables in the inverse model, {\eqn{n x p}} matrix, each row is a response vector
}
\item{Fy}{vector of covariates in the inverse problem, vector containing functions of the response variable in the original multiple regression problem. Is a {\eqn{n x p}} matrix, each row is the corresponding response vector
}
  \item{aux}{list containing the constants for the tau-scale (with components c1, k1, c2, k2, as described in  %\code{\link[kappa_and_c]}{kappa_and_c}
  \code{\link[tauPFC]{kappa_and_c}})
  %\href{kappa_and_c}{kappa_and_c})
}
  \item{efficiency}{required efficiency for initial robust estimator
}
}
\details{This function computes the coefficient matrix concatenating the coefficients resulting from the fitting of univariate regressions using a robust estimator with high breakdown point and high efficiency (require "robust"). Univariate estimators are computed using  %\href{robust::lmRob.Rd}{lmRob}
\code{\link[robust]{lmRob}}
"lmRob" function. It computes a robust estimator of errors covariance matrix (require "rrcov").
}
\value{
The estimation is given as a list with components
\item{beta0 }{coefficient matrix }
\item{delta0 }{S covariance matrix of residuals}
}
\seealso{For the final estimation use \code{\link{tauestimate}} with this initial estimator as an input
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
