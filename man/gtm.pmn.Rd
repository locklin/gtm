\name{gtm.pmn}
\alias{gtm.pmn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculates the posterior mean projection of data into the latent space }
\description{
 The posterior mean projection of a point from the target space, t, is the mean 
 of the corresponding posterior distribution induced in the latent space. 
}
\usage{
gtm.pmn(TT, X, FI, W, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TT}{ data points representing the distribution in 
     the target space. N-by-D  }
  \item{X}{ data points forming a latent variable sample 
     of the distribution in the latent space. K-by-L}
  \item{FI}{ activations of the basis functions when fed X; K-by-(M+1)  }
  \item{W}{ a matrix of trained weights }
  \item{b}{ the trained value for beta}
}
\value{
\code{output} - the posterior means in latent space. N-by-L 
}
\seealso{\code{\link{gtm.ppd1}},\code{\link{gtm.ppd2}} \code{\link{gtm.pmd}}}
\examples{
   lv <- gtm.demo()
   gtm.pmn(lv$TT, lv$X, lv$FI, lv$W, lv$beta)
 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{models}
