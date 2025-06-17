\name{easy}
\alias{easy}
\title{An example with M=5 (basis functions) and r=3 (non-zero eigenvalues)}
\description{A simulated dataset as an example which corresponds to the "easy" case in the paper}
\format{
easy is a list with six components (in the given order):
\describe{
\item{data}{data matrix with three columns: column 1--ID, column 2--measurement, column 3--time.}
\item{eigenfunctions}{true eigenfunctions: generated from cubic Bsplines with M=5 equally spaced knots.}
\item{eigenvalues}{true eigenvalues: first--1, second--0.66, third--0.52, others--zero.}
\item{number_of_basis}{true number of basis functions: M=5.}
\item{dimension }{true dimension of the process: r=3.}
\item{error_sd }{true error standard deviation: 0.25.}
 }
}

\details{mean curve of the process is zero; principal component scores and errors are all i.i.d N(0,1); there are 200 subjects, and each has 
2~10 measurements uniformly distributed on [0,1]; in total there are 1227 measurements}
\usage{data(easy)}
\keyword{datasets}