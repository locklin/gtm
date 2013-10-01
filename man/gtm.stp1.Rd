\name{gtm.stp1}
\alias{gtm.stp1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generates the components of a GTM with a 1D latent space}
\description{
 1-5 lines 
}
\usage{
gtm.stp1(T, noLatVarSmpl, noBasFn, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{target data, to be modelled by the GTM.}
  \item{noLatVarSmpl}{number of samples in the latent variable space}
  \item{noBasFn}{number of basis functions}
  \item{s}{the width of basis functions relative
				to the distance between two neighbouring 
				basis function centres, i.e. if s = 1,
				the basis functions will have widths
				(std.dev) equal to (i.e. 1 times) the  
				distance between two neighbouring basis  
				function centres. }
}
\value{
  List with the following components
  \item{\code{X}}{the grid of data points making up the latent variable sample;
      a vector of length noLatVarSmpl, in which each row is a data point}
  \item{\code{MU}}{a noBasFn-element vector holding the coordinates of the centres
      of the basis functions }
  \item{\code{FI}}{the activations of the basis functions when fed the latent variable
      sample X, and a bias unit fixed to 1.0; a matrix with the same number 
      of rows as X and noBasFn+1 columns (+1 for the bias)  }
  \item{\code{W}}{ the initial matrix of weights, mapping the latent variable sample X
      linearly onto the first principal component of the target data (T) }
  \item{\code{beta}}{- the intial value for the inverse variance of the data space
      noise model  }
}
\note{ 
The latent variable sample is constructed as a uniform grid on the 
interval [-1, 1]. Similarly the centres of the basis function are
 gridded uniformly over the latent variable sample, with equal 
 standard deviation, set relative to the distance between two centres. 
The initial linear mapping maps mean and std.devs.
 1:1 from the latent to the target sample along the principal component. 
}
\seealso{ \code{\link{gtm.stp2}} }
\examples{
  T = matrix(3:61/20, ncol=1);
  T = cbind(T, T + 1.25 * sin(2*T));
  # plot of sample data 
  plot(T[,1], T[,2], col='red', pch=21, xlim=c(0, 3.5), ylim=c(0,3.5))
    
  # setup and training
  stpResult = gtm.stp1(T, 20, 5, 2);
  stpY = stpResult$FI %*% stpResult$W
  
  # plot initialized data points
  matplot(stpY[,1], stpY[,2], col='green', add=TRUE, type='l')
  matplot(stpY[,1], stpY[,2], col='green', add=TRUE, type='p', pch=23)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}
