\name{gtm.pci}
\alias{gtm.pci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Returns a weight matrix initialised using principal components. }
\description{
 The returned weight matrix maps the mean of the latent variable to the mean of
  the target variable, and the L-dimensional latent variable variance to the 
  variance of the target data along its L first principal components.

}
\usage{
W = gtm.pci(TT, X, FI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TT}{target distribution sample; one data point per row; N-by-D  }
  \item{X}{ the latent distribution sample, K-by-L  }
  \item{FI}{basis functions' activation when fed the latent data,
     X, plus a bias, K-by-(M+1)  }
}
\value{
 \code{W} - the initialised weight matrix, (K+1)-by-D 
}

\note{
   The first dimension (column) of X will map to the first principal
 component, the second dimension (column) of X will map to the 
 second principal component, and so on. This may be of importance 
 for the choice of sampling density along the different dimensions
  of X, if it differs between different dimensions.
}
\seealso{ \code{\link{gtm.ri}}, \code{\link{gtm.pci.beta}} }
\examples{
  MU = matrix(c(0,1), ncol = 1)
  sigma = 1
  X = gtm.pts(8)
  FI = gtm.gbf(MU, sigma, X)
  t = matrix((1:40)/5, ncol=1)
  TT = cbind(t, sin(t))
  gtm.pci(TT, X, FI)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}
