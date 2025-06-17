\name{fpca.score}
\alias{fpca.score}
\title{Functional principal component scores}
\description{A function to estimate the functional principal component scores by the best linear unbiased predictors (Yao et al. 2005). 
}
\usage{
fpca.score(data.m,grids.u,muhat,eigenvals,eigenfuncs,sig2hat,K)
}

\arguments{
    \item{data.m}{Matrix with three columns. Each row corresponds to one measurement for one subject.
     Column One: subject ID (numeric or string); 
     Column 2: measurement (numeric); Column 3: corresponding measurement time (numeric);
     Missing values are not allowed. Same format as the data input for \link{fpca.mle}.
     }
     
     \item{grids.u}{Grid of time points used in evaluating the mean and eigenfunctions (on the original scale). Same as 'grid' returned by \link{fpca.mle}.
     }
     
     \item{muhat}{Mean evaluated on the same grids as in grids.u. An estimate is returned by \link{fpca.mle}.
     }
     
     \item{eigenvals}{Eigenvalues. An estimate is returned by \link{fpca.mle}.
     }
      
     \item{eigenfuncs}{Eigenfunctions evaluated on the same grids as in grids.u. An estimate is returned by \link{fpca.mle}.
     }
      
     \item{sig2hat}{Noise variance. An estimate is returned by \link{fpca.mle}.
     }
     
     \item{K}{Number of eigenfunctions used to derive the fpc scores.
     }
     
  
   } 
\details{
    'fpca.score' uses best linear unbiased predictors (BLUP) to estimate the functional principal component scores for each subject
    } 

\value{
    An n by K matrix containing the first K functional principal component scores for each subject.  
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
\link{fpca.mle} for model fitting, \link{fpca.pred} for prediction
}

\author{J. Peng, D. Paul}

\keyword{methods}


