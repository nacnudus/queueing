% File man/VNq.o_MMCK.Rd
\name{VNq.o_MMCK}
\alias{VNq.o_MMCK}
\title{Returns the variance of the number of customers in the queue in the M/M/c/K queueing model}
\description{
  Returns the variance of the number of customers in the queue in the M/M/c/K queueing model
}
\usage{
  \method{VNq}{o_MMCK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCK}
  \item{\dots}{aditional arguments}
}

\details{Returns the variance of the number of customers in the queue in the M/M/c/K queueing model}

\references{
  [Sztrik2012] Dr. Janos Sztrik (2012).\cr
\emph{Basic Queueing Theory}.\cr
University of Debrecen, Faculty of Informatics.
}

\seealso{
  \code{\link{QueueingModel.i_MMCK}}.
}


\examples{

## create input parameters
i_mmck <- NewInput.MMCK(lambda=8, mu=4, c=5, k=12)

## Build the model
o_mmck <- QueueingModel(i_mmck)

## Returns the variance
VNq(o_mmck)
}

\keyword{M/M/c/K}
