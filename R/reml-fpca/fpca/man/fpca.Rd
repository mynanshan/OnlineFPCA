\name{fpca}
\alias{fpca}
\title{The fpca package: summary information}
\description{The package implements the restricted maximum likelihood estimation through a Newton-Raphson procedure described in Peng and Paul (2009) 
to estimate functional principal components from (sparsely and irregularly observed) longitudinal data.
}
\details{This is version 0.2-1 updated in Feb, 2011.  Two new functions, fpca.score, fpca.pred, are included. 
Missing values are not allowed. Subjects with only one measurement will be automatically excluded. 
The main function is 'fpca.mle'. Simulated data sets can be called by 'data(easy)' and 'data(prac)'. Type 'help(easy)' and 'help(prac)' 
to see details.
Packages 'sm' and 'splines' are used by this package.  The code for EM (as initial estimate) is provided by Professor G. James in USC (with slight modifications).
}

\references{
Peng, J. and Paul, D. (2009).  
A geometric approach to maximum likelihood estimation of the functional principal components from sparse longitudinal data.
 Journal of Computational and Graphical Statistics. December 1, 2009, 18(4): 995-1015 

James, G. M., Hastie, T. J. and Sugar, C. A. (2000) Principal component models for sparse
functional data. Biometrika, 87, 587-602.

Yao, F., Mueller, H.-G. and Wang, J.-L. (2005) Functional data analysis for sparse longitudinal
data. Journal of the American Statistical Association 100, 577-590

}
\seealso{
\link{fpca.mle} for model fitting, \link{fpca.score} for fpc scores, \link{fpca.pred} for prediction
}

\author{J. Peng, D. Paul}
\keyword{package}