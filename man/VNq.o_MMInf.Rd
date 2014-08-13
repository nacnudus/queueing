% File man/VNq.o_MMInf.Rd
\name{VNq.o_MMInf}
\alias{VNq.o_MMInf}
\title{Returns the variance of the number of customers in the queue in the M/M/Infinite queueing model}
\description{
  Returns the variance of the number of customers in the queue in the M/M/Infinite queueing model
}
\usage{
  \method{VNq}{o_MMInf}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Returns the variance of the number of customers in the queue in the M/M/Infinite queueing model}

\references{
 [Sztrik2012] Dr. Janos Sztrik (2012).\cr
\emph{Basic Queueing Theory}.\cr
University of Debrecen, Faculty of Informatics.
}

\seealso{
  \code{\link{QueueingModel.i_MMInf}}.
}


\examples{
## create input parameters
i_mminf <- NewInput.MMInf(lambda=0.25, mu=4, n=0)

## Build the model
o_mminf <- QueueingModel(i_mminf)

## Returns the variance
VNq(o_mminf)
}

\keyword{M/M/Infinite}
