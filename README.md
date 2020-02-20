# tauPFC
Computes robust estimators for the PFC model

This repo holds a package in R that compute robust estimators for the sufficient dimension
reduction problem. It estimates the parameters of 
the Principal Fitted Components (PFC) model, using tau estimators that can resist the presence of outliers. 
We focus on the reduction subspace estimation. As PFC model and multivariate reduced-rank regression are closely
related, our proposal can also be used in the estimation of this model.
The package also hold a routine to perform classical estimation (MLE), 
as proposed by Cook and Forzani (2008).
We further propose a cross validation
procedure to select the dimension of the subspace.

The complete description and theoretical
aspects of the robust estimation technique proposed can be found in Bergesio, A., Szretter Noste, M.E. and Yohai, V.J.,
"A robust proposal of estimation for the sufficient dimension reduction problem".


# Installation

You can install pkgreviewr from GitHub with:
```
# install.packages("devtools") 
devtools::install_github("meszre/tauPFC")
```
