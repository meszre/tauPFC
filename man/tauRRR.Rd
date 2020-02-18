\name{tauRRR}
\alias{tauRRR}
\title{Tau estimate for Reduced Rank Regression
}
\description{Fits a reduced-rank regression (RRR) model (for a description, see Izenman(2008)) with a robust procedure that can resist the presence of outliers. It computes tau type estimates, see Bergesio et al. (2020). It fits also a multivariate linear model (MLM) without rank restriction if the rank is chosen properly.
}
\usage{
tauRRR(yy, XX, d, efficiency = 0.90)
}
\arguments{
  \item{yy}{vector of response variables, {\eqn{n x p}} matrix, each row is a response vector
}
\item{XX}{vector of covariates in the multivariate regression problem. Is a {\eqn{n x r}} matrix, each row is the corresponding covariate vector}


\item{d}{rank of coefficient matrix. If \code{d}\eqn{ = min(p,r)}, this routine produces robust multivariate linear regression estimators. If \code{d}\eqn{ < min(p,r)}, this routine produces a reduced rank regression model, with arbitrary covariance of errors }

\item{efficiency}{efficiency of the robust estimators, 0.90 by default}
}


\details{We consider the multivariate linear reduced rank regression (RRR) model given by
\deqn{Y = \mu + ABX + \epsilon,}
where the variables are
\itemize{
\item{} \eqn{Y} is a \emph{p x 1} observed response vector
\item{} \eqn{X} is a \emph{r x 1} observed covariates vector
\item \eqn{\epsilon} is unobserved  \emph{p} dimensional vector of errors}
and the unknown parameters (to be estimated) are
\itemize{
\item \eqn{\mu} a \emph{p x 1} vector of intercepts
\item{} \emph{A} is a full-rank \emph{p x d} matrix
\item{} \emph{B} is a full-rank \emph{d x r} matrix
\item{} cov({\eqn{\epsilon}}), the \emph{p x p} covariance matrix of errors {\eqn{\epsilon}}, }
Both coefficient matrices \emph{A} and \emph{B} are not unique, but their product \emph{p x r} matrix \emph{C = AB} is unique. The rank of \emph{C}  is d \eqn{\le min(p,r)}. The notation refers to Izenman (2008).}
\value{List with the following components, of the RRR model described above
\item{mu}{tau-estimator for intercept vector}
\item{RRRcoef}{tau-estimator for the coefficient \emph{p x r} matrix \emph{C} of rank \code{d} }
\item{AA}{tau-estimator for \emph{A}, a full-rank \emph{p x d} matrix}
\item{BB}{tau-estimator for \emph{B}, a full-rank \emph{d x r} matrix}
\item{cov.error}{tau-estimator for the covariance matrix of errors, cov({\eqn{\epsilon}})}
}

\references{
Izenman, A. J. (2008). Modern multivariate statistical techniques. {\emph{Regression, classification and manifold learning}}, New York: Springer.

Bergesio, A., Szretter Noste, M. E. and Yohai, V. J. (2020). A robust proposal of estimation for the
sufficient dimension reduction problem
}

\examples{

# We work with Example 6.3.3, Chemical Composition of Tobacco, from
# Izenman (2008). Dataset available in \link[rrr]{tobacco} or in
# \url{https://astro.temple.edu/~alan/tobacco.txt}

library(rrr)
data(tobacco)

XX = as.matrix(tobacco[,4:9])  # covariates
yy = as.matrix(tobacco[,1:3])  # response vector

###############################
# RRR model with d=2
###############################
# robust MLM fit
robustRRR2 = tauRRR(yy, XX, d=2)

# maximum likelihood MLM, with arbitrary covariance of errors matrix
maxliRRR2 = MLE(yy, XX, d=2)


###############################
# robust MLM fit
###############################
robustMLM = tauRRR(yy, XX, d=3)

# classical MLM, with covariance of errors a multiple of identity
classicalMLM = lm(yy ~ XX)

# maximum likelihood MLM, with arbitrary covariance of errors matrix
maxlikMLM = MLE(yy, XX, d=3)

# to show that the three of them agree, we can do any of the following

# 1. Verify they span the same column space
library(pracma)
angle_AB(orth(robustMLM$RRRcoef),orth(t(classicalMLM$coefficients[2:7,])))
# equivalently
angle_AB(robustMLM$AA ,orth(t(classicalMLM$coefficients[2:7,])))
angle_AB(robustMLM$AA,maxlikMLM$gamma)

# 2. Plot coefficients estimated by every method
plot(robustMLM$RRRcoef,t(classicalMLM$coefficients[2:7,]))
abline(0,1)

points(t(classicalMLM$coefficients[2:7,]),maxlikMLM$gamma\%*\%maxlikMLM$beta,col="red")

# 3. Compute maximum absolute difference in coefficients estimated by every method
max(abs(t(classicalMLM$coefficients[2:7,])-maxlikMLM$gamma\%*\%maxlikMLM$beta))
max(abs(t(classicalMLM$coefficients[2:7,])-robustMLM$RRRcoef))
}

