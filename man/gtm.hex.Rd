\name{gtm.hex}
\alias{gtm.hex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produces a 2D grid with points arranged in a hexagonal lattice. }
\description{
  The grid is centered on the origin and scaled so the dimension (X or Y) 
   with largest number of points ranges from -1 to 1. 
}
\usage{
gtm.hex(xdim, ydim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdim, ydim}{number of points along the X and Y dimensions, 
     respectively; must be >=2.  }
 
}
\value{
  \code{grid} - a (xdim*ydim)-by-2 matrix of grid points with the first point being 
  the top-left corner and subsequent points following columnwise. 
}
\seealso{\code{\link{gtm.rctg}} }
\examples{
  plot(gtm.hex(3,7))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}