\name{gtm.lbf}
\alias{gtm.lbf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculates the output of linear basis functions for a given set of inputs }
\description{
  This simply amounts to returning the set of inputs with an extra bias column
   of ones after the last column in the input set matrix. 
}
\usage{
gtm.lbf(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ the latent variable sample forming the set of inputs; K-by-L }
}

\value{
  \code{FI} - the matrix of basis functions output values; K-by-(L+1),
   "+1" for the bias basis function 
}

\seealso{ \code{\link{gtm.gbf}} }
\examples{
  gtm.lbf(matrix(c(3,4,2,7), ncol=2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}