\name{kappa_and_c}
\alias{kappa_and_c}

\title{Constants for tau-scale estimator}
\description{This function computes the constants to define the rho-functions (in the rho quasi-optimal family) and expectations in order to obtain high efficiency and the maximum possible breakdown point (0.5).
The values are obtained by solving the equations that relate the rho-functions with these characteristics of robustness}
\usage{kappa_and_c(p, efficiency)}
\arguments{
  \item{p}{dimension of the vector of the response variable}
  \item{efficiency}{efficiency required, values between 0 and 1, excluding extremes}}
\details{
For more details about the way the constants are chosen, see Bergesio et al. (2020). For the definition of rho quasi optimal family of functions, see Muler and Yohai (2002) and Salibian-Barrera et al. (2016).
}
\value{
The result is returned as a list containing
 \item{c1}{ tuning parameter such that {\eqn{\rho}}{\out{<sub>1</sub>}}  is strictly increasing in {\out{[0; c<sub>1</sub>]}} and constant for values greater than {\out{ c<sub>1</sub>}}.}
  \item{c2}{ tuning parameter such that {\eqn{\rho}}{\out{<sub>2</sub>}} is strictly increasing in {\out{[0; c<sub>2</sub>]}} and constant for values greater than {\out{ c<sub>2</sub>}}.}
 \item{k1}{ constant value involved in the {\eqn{\rho}}{\out{<sub>1</sub>}} expectation.}
 \item{k2}{ constant value involved in the {\eqn{\rho}}{\out{<sub>2</sub>}} expectation.}
 }

\references{
Bergesio, A., Szretter Noste, M. E. and Yohai, V. J. (2020). A robust proposal of estimation for the
sufficient dimension reduction problem.

Muler,  N. and Yohai,  V. J. (2002). Robust estimates for arch processes.
{\emph{Journal of Time Series Analysis}}, 23(3):341-375.

Salibian-Barrera, M., Van Aelst, S. and Yohai, V. J. (2016). Robust tests for linear regression models based on \eqn{\tau}-estimates. {\emph{Computational Statistics & Data Analysis}} 93:436-455.
}



\note{It is not an efficient procedure, so it can take some seconds.}

\examples{
kappa_and_c(p=1, efficiency=0.95)}
