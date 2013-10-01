\name{gtm.trn}
\alias{gtm.trn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Optimize (train) the parameters of a GTM model, using an EM algorithm.}
\description{
The function first computes distances of points in data space and corresponding
conditional probabilities (called responsibilities). It then uses analytic
formulas to compute new values of W and beta. The formulas involve inverting
a matrix, which is carried out by QR algorithm, if the matrix is positive 
definite, or by computing its pseudo-inverse, if the matrix has zero eigenvalue.
}
\usage{
gtm.trn(T, FI, W, l, cycles, beta, m = 1, quiet = FALSE, minSing = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{matrix containing a sample of the distribution 
		to be modeled; N-by-D}
  \item{FI}{matrix containing the output values from
			the basis functions, when fed the latent
			variable sample; K-by-(M+1)}
  \item{W}{ an initial weight matrix; (M+1)-by-D}
  \item{l}{ weight regularisation factor}
  \item{cycles}{ number of training cycles} 
  \item{beta}{an initial value for beta, the inverse variance of 
			the Gaussian mixture generated in the data space}
  \item{m}{ mode of calculation; it can be set to 0, 1 or 2 
			corresponding to increasingly elaborate 
			measure taken to reduce the amount of
			numerical errors; mode = 0 will be fast but 
			less accurate, mode = 2 will be slow but 
			more accurate; the default mode is 1 }
  \item{quiet}{quiet execution; if q equals to TRUE,
			the plotting and echoing of the values of log- 
			likelihood and beta during training is supressed. 
			This argument is optional; if omitted the training 
			is run non-quiet. }
  \item{minSing}{ cutoff parameter controlling calculation of the pseudo-inverse matrix. Singular values less than
\code{minSing} are considered numerically zero. }
}
\value{
List with the following components
\item{\code{W}, \code{beta}} the corresponding weight matrix and inverse variance after training
\item{llhLog} the log-likelihood after each cycle of training; optional output argument
}

\examples{
  T = matrix(3:61/20, ncol=1);
  T = cbind(T, T + 1.25 * sin(2*T));
  # plot of sample data 
  plot(T[,1], T[,2], col='red', pch=21, xlim=c(0, 3.5), ylim=c(0,3.5))
    
  # setup and training
  stpResult = gtm.stp1(T, 20, 5, 2);
  trnResult = gtm.trn(T, stpResult$FI, stpResult$W, 0.0, 20, stpResult$beta);

  # plot of trained results
  Y = stpResult$FI %*% trnResult$W;
  plot(Y[,1], Y[,2], col='green', type='l')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}
