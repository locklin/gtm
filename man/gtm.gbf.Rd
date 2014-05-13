\name{gtm.gbf}
\alias{gtm.gbf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculates the output of Gaussian basis functions for a given set of input  }
\description{First the function computes matrix \eqn{(d.{ij})} of  
distances between centers  of
Gaussians given by \code{MU} and latent points \code{X}. The output values are 
given  by \deqn{\Phi.{ij} = \exp(- \frac{d.{ij}}{2 \sigma^2})}{%
             p(x) = exp(-1/(2* sigma^2)*d.{ij})}
	     }

\usage{
gtm.gbf(MU, sigma, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MU}{ a M-by-L matrix containing the centers of the basis functions }
  \item{sigma}{ a scalar giving the standard deviation of the radii-symmetric Gaussian basis functions }
  \item{X}{ the latent variable sample forming the set of inputs; K-by-L  }
}

\value{
 \code{output} - the matrix of basis functions output values; K-by-(M+1), "+1" for a bias basis function with a fixed value of 1.0 
}
\seealso{  \code{\link{gtm.lbf}} }
\examples{  
  MU <- matrix(c(0,1), ncol = 1)
  sigma <- 1
  X <- gtm.pts(8)
  FI <- gtm.gbf(MU, sigma, X)
  matplot(X, FI, type="l")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

