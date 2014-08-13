% File man/VNq.o_MMCKM.Rd
\name{VNq.o_MMCKM}
\alias{VNq.o_MMCKM}
\title{Returns the variance of the number of customers in the queue in the M/M/c/K/m queueing model}
\description{
  Returns the variance of the number of customers in the queue in the M/M/c/K/m queueing model
}
\usage{
  \method{VNq}{o_MMCKM}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCKM}
  \item{\dots}{aditional arguments}
}

\details{Returns the variance of the number of customers in the queue in the M/M/c/K/m queueing model}

\references{
  [Sztrik2012] Dr. Janos Sztrik (2012).\cr
\emph{Basic Queueing Theory}.\cr
University of Debrecen, Faculty of Informatics.
}

\seealso{
  \code{\link{QueueingModel.i_MMCKM}}.
}


\examples{
## create input parameters
i_mmckm <- NewInput.MMCKM(lambda=0.25, mu=4, c=2, k=4, m=8, method=0)

## Build the model
o_mmckm <- QueueingModel(i_mmckm)

## Returns the variance
VNq(o_mmckm)
}

\keyword{M/M/c/K/m}
