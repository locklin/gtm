\name{gtm.resp6}
\alias{gtm.resp6}
\alias{gtm.resp3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Log-likelihood and component responsibilities under a Gaussian mixture  }
\description{
The responsibility R(k,n) is the probability of a particular component in the
 Gaussian mixture, k, having generated a particular data point, n. It is 
 calculated from the distances between the data point n and the centres of 
 the mixture components, 1..K, and the inverse variance, beta, common to 
 all components. 
}
\usage{
gtm.resp6(DIST, minDist, maxDist, beta, D, mode)
gtm.resp3(DIST, beta, D)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DIST}{a K-by-N matrix in which element (k,n) is the squared distance 
      between the centre of component k and the data point n.  }
  \item{minDist, maxDist}{ vectors containing the minimum and maximum of 
      each column in DIST, respectively; 1-by-N; required iff m > 0}
  \item{beta}{ a scalar value of the inverse variance common to all components 
      of the mixture}
  \item{D}{ dimensionality of space where the data and the Gaussian mixture 
      lives; necessary to calculate the correct log-likelihood }
  \item{mode}{optional argument used to control the mode of calculation; it can
     be set to 0, 1 or 2 corresponding to increasingly elaborate measure taken 
     to reduce the amount of numerical errors; mode = 0 will be fast but less
      accurate, mode = 2 will be slow but more accurate; the default mode is 0}
}
\value{
  A list containing 
  \item{\code{llh}}{the log-likelihood of data under the Gaussian mixture }
  \item{\code{R} }{an K-by-N responsibility matrix; R(k,n) is the responsibility
     takened by mixture component k for data point n. }

}
\seealso{ \code{\link{gtm.dist}} }
\examples{
  lv = gtm.demo()
  gtm.resp3(gtm.dist(lv$T, lv$Y), lv$beta, 2)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}
