\name{gtm.rctg}
\alias{gtm.rctg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produces a 2D grid with points arranged in a rectangular lattice }
\description{
  The grid is centered on the origin and scaled so the dimension (X or Y) 
  with largest number of points ranges from -1 to 1. 
}
\usage{
grid = gtm.rctg(xDim, yDim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xDim,yDim}{ number of points along the X and Y dimensions, 
       respectively; must be >=2.  }
}
\value{
\code{grid} - a (xDim*yDim)-by-2 matrix of grid points with the first point 
  being the top-left corner and subsequent points following column-wise. 
}

\seealso{\code{\link{gtm.hxg}}}
\examples{
  plot(gtm.rctg(3,4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}
