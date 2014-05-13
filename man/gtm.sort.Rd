\name{gtm.sort}
\alias{gtm.sort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sorts the columns of argument matrix R in increasing order }
\description{
Size of the column is the usual Euclidean norm.
}
\usage{
gtm.sort(R)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{ an (unsorted) matrix  }
}

\value{
\code{output} - the corresponding sorted matrix
}
\examples{
  m <- matrix( c(1,4,2,7,1,3), nrow=2)

  gtm.sort(m)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}
