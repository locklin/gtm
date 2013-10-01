\name{gtm-internal}
\alias{gtm-internal}
\alias{gtm.r2m1}
\alias{gtm.r2m2}
\alias{gtm.r2m3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Internal functions in gtm package}
\description{
Utility functions not intended for common usage.
}
\usage{
gtm.r2m1 (cX, meshRows, meshCols)
gtm.r2m2 (cX, cY, meshRows, meshCols)
gtm.r2m3 (cX, cY, cZ, meshRows, meshCols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{cX, cY, cZ}{ column vectors with x-, y-, and x-data respectively; N-by-1}

\item{meshRows, meshCols}{number of rows and colmuns of
	 the mesh matrices; meshRows*meshCols = N }
}

\value{
List of mesh matrices.
}
\keyword{models}

