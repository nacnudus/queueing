% File man/VT.o_MMC.Rd
\name{VT.o_MMC}
\alias{VT.o_MMC}
\title{Returns the variance of the time spend in the M/M/c queueing model}
\description{
  Returns the variance of the time spend in the M/M/c queueing model
}
\usage{
  \method{VT}{o_MMC}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMC}
  \item{\dots}{aditional arguments}
}

\details{Returns the variance of the time spend in the M/M/c queueing model}

\references{
  [Sztrik2012] Dr. Janos Sztrik (2012).\cr
\emph{Basic Queueing Theory}.\cr
University of Debrecen, Faculty of Informatics.
}

\seealso{
  \code{\link{QueueingModel.i_MMC}}.
}


\examples{

## create input parameters
i_mmc <- NewInput.MMC(lambda=5, mu=10, c=2, n=0, method=0)

## Build the model
o_mmc <- QueueingModel(i_mmc)

## Returns the variance of the time spend in the system
VT(o_mmc)
}

\keyword{M/M/c}
