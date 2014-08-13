% File man/L.o_MMInfKK.Rd
\name{L.o_MMInfKK}
\alias{L.o_MMInfKK}
\title{Returns the mean number of customers in the M/M/Infinite/K/K queueing model}
\description{
  Returns the mean number of customers in the M/M/Infinite/K/K queueing model
}
\usage{
  \method{L}{o_MMInfKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMInfKK}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean number of customers in the M/M/Infinite/K/K queueing model}

\references{
  [Kleinrock1975] Leonard Kleinrock (1975).\cr
\emph{Queueing Systems Vol 1: Theory}.\cr
John Wiley & Sons.
}

\seealso{
  \code{\link{QueueingModel.i_MMInfKK}}.
}


\examples{
## create input parameters
i_MMInfKK <- NewInput.MMInfKK(lambda=0.25, mu=4, k=4)

## Build the model
o_MMInfKK <- QueueingModel(i_MMInfKK)

## Returns the L
L(o_MMInfKK)
}

\keyword{M/M/Infinite/K/K}
