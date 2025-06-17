\name{fpca.pred}
\alias{fpca.pred}
\title{Predicted trajectories}
\description{A function to predict trajectory for each subject 
}
\usage{
fpca.pred(fpcs, muhat,eigenfuncs)
}

\arguments{
    \item{fpcs}{Functional principal component scores. An estimate is returned by \link{fpca.score}}
    \item{muhat}{Mean curve evaluated on a grid. An estimate is returned by \link{fpca.mle}.
     }
     
    \item{eigenfuncs}{Eigenfunctions evaluated on the same grids as in 'muhat'. An estimate is returned by \link{fpca.mle}.
     } 
     
   } 
\details{
    'fpca.pred' gives predicted trajectories (evaluated on a fine grid). 
    } 

\value{
    A matrix where each column corresponds to the predicted trajectory of a subject.
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
\link{fpca.mle} for model fitting, \link{fpca.score} for fpc scores
}

\author{J. Peng, D. Paul}

\keyword{methods}


