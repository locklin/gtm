\name{gtm.ri}
\alias{gtm.ri}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns an initial random weight matrix }
\description{
  Generates a weight matrix with the bias weights set to map to the mean 
  of the target distribution and remaining weights drawn at random from 
  a Gaussian distribution with zero mean and variances choosen so that 
  the variances of the generated distribution roughly match the variances 
  of the target distribution.

  In addition, an initial value for beta may be calculated as the inverse 
  of the average distance between each Gaussian centre, calculated with the 
  random mapping, and its nearest neighbours in the set of data points. 
}
\usage{
gtm.ri(TT, FI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TT}{sample of target distribution, used for calculating the mean; 
     one data point per row; N-by-D }
  \item{FI}{ basis functions' activations when fed the latent data, X, 
       plus a bias, K-by-(M+1)  }
}

\value{
 \code{output} - the initialised weight matrix; (M+1)-by-D 
}

\seealso{\code{\link{gtm.pci}},\code{\link{gtm.pci.beta}} }
\examples{
  TT <- matrix(3:61/20, ncol=1);
  TT <- cbind(TT, TT + 1.25 * sin(2*TT));
  # plot of sample data 
  plot(TT[,1], TT[,2], col='red', pch=21, xlim=c(0, 3.5), ylim=c(0,3.5))
    
  # setup and training
  X <- gtm.pts(20)
  MU <- gtm.pts(5)
  sigma <- 1
  FI <- gtm.gbf(MU, sigma, X)
  W <- gtm.ri(TT, FI);
  Y <- FI %*% W;
  
  # plot of initialized data
  matplot(Y[,1], Y[,2], col='green', add=TRUE, type='p', pch=23)
  beta <- gtm.bi(Y)
  symbols(x=Y[,1], y=Y[,2], circles=matrix(sqrt(1/beta) ,ncol=1,nrow=nrow(Y)),
     fg='green',add=TRUE, inches=FALSE);
  
  trnResult <- gtm.trn(TT, FI, W, 0.0, 20, beta, quiet=TRUE);

  # plot of trained results
  Y <- FI %*% trnResult$W;
  matplot(Y[,1], Y[,2], col='blue', type='l', add=TRUE)
  matplot(Y[,1], Y[,2], col='blue', type='p', add=TRUE, pch=22)
  b <- trnResult$beta
  symbols(x=Y[,1], y=Y[,2], circles=matrix(sqrt(1/b) ,ncol=1,nrow=nrow(Y)),
     fg='blue',add=TRUE, inches=FALSE);
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}
