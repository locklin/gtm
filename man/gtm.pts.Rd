\name{gtm.pts}
\alias{gtm.pts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate division of interval}
\description{
Returns M points equally spaced in the interval [-1,1].
}
\usage{
gtm.pts(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{ the number of points }
}

\value{
  Returns M-by-1 matrix. In the first row is -1, in the last row 1.  
}

\examples{
  gtm.pts(9)
}
\keyword{cluster}
\keyword{models}
