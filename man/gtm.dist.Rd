\name{gtm.dist}
\alias{gtm.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate the squared distances between two sets of data points}
\description{
  This function calculates distances between all data points in the two 
    data sets TT and Y and returns them in a matrix. 
}
\usage{
gtm.dist(x, y, m = 0, kind="Euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y}{data set matrices in which each row is a data point; 
     dimensions N-by-D and K-by-D respectively }
  \item{m}{mode of calculation; iff m > 0, min- and maxDist (below) are
    calculated; the default mode is0  }
  \item{kind}{ type of distance; presently only tested with Euclidean,
    but could eventually include all of the proxy methods in
    gtm::GLOBALPROXYNAMES }

}
\details{
 
}
\value{
if m == 0 then \code{DIST} - matrix containing the calculated distances; 
    dimension K-by-N; DIST(k,n) contains the squared distance between TT(n,:) and Y(k,:). 
 if m > 0 then a list containing matrix \code{DIST}, \code{minDist}, 
  \code{maxDist}, the latter two being 
  vectors containing the minimum and maximum of each column in DIST, respectively; 1-by-N
}
\examples{  
  A <- matrix(c(1,3,4,2,0,-1), ncol=2)
  B <- matrix(c(1,0,0,1), ncol=2)
  gtm.dist(A, B)
  gtm.dist(A, B ,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}
