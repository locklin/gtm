\name{gtm.demo}
\alias{gtm.demo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Demonstrates the GTM with a 2D target space and a 1D latent space}
\description{
This script generates a simple data set in 2 dimensions, with an intrinsic dimensionality of 1, and trains a GTM with a 1-dimensional latent variable to model this data set, visually illustrating the training process. 
}
\usage{
gtm.demo(randomSetup = FALSE, points = 5, samples = 20,interactive=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{randomSetup}{if FALSE, use \code{\link{gtm.pci.beta}} for initial configuration, otherwise, initialize using \code{\link{gtm.ri}}  }
  \item{points}{number of Gaussians}
  \item{samples}{number of latent space samples}
  \item{interactive}{set to TRUE by default; turn off for timing/benchmarks}
}

\value{
 List with elements
  \item{\code{X}}{latent sample points}
  \item{\code{TT}}{data in the target space}
  \item{\code{W}}{the initialised weight matrix}
  \item{\code{Y}}{centers of Gaussians}
  \item{\code{beta} }{the optimized beta value, scalar. }
  \item{\code{llh}}{log-likelihood}
  \item{\code{FI}}{basis functions' activation}
}
\note{
}
\seealso{  \code{\link{gtm.ri}}, \code{\link{gtm.pci}} }
\examples{
  gtm.demo ()
  gtm.demo (randomSetup = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}
