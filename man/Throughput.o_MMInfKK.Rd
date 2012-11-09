% File man/Throughput.o_MMInfKK.Rd
\name{Throughput.o_MMInfKK}
\alias{Throughput.o_MMInfKK}
\title{Throughput of a M/M/Infinite/K/K queueing model}
\description{
  Returns the throughput of a M/M/Infinite/K/K queueing model
}
\usage{
  \method{Throughput}{o_MMInfKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMInfKK}
  \item{\dots}{aditional arguments}
}

\details{Returns the throughput of a M/M/Infinite/K/K queueing model}

\references{
  [Kleinrock1975] Leonard Kleinrock (1975).\cr
\emph{Queueing Systems Vol 1: Theory}.\cr
John Wiley & Sons.
}

\seealso{
  \code{\link{NewInput.MMInfKK}}, \code{\link{CheckInput.i_MMInfKK}}, \code{\link{QueueingModel.i_MMInfKK}}
}


\examples{
## create input parameters
i_MMInfKK <- NewInput.MMInfKK(lambda=0.25, mu=4, k=4)

## Build the model
o_MMInfKK <- QueueingModel(i_MMInfKK)

## Throughput
Throughput(o_MMInfKK)

}
\keyword{M/M/Infinite/K/K}
