% File man/summary.o_MMInfKK.Rd
\name{summary.o_MMInfKK}
\alias{summary.o_MMInfKK}
\title{Reports the results of a M/M/Inf/K/K queueing model}
\description{
  Reports the results of a M/M/Inf/K/K queueing model.
}
\usage{
  \method{summary}{o_MMInfKK}(object, \dots)
}
\arguments{
  \item{object}{a object of class o_MMInfKK}
  \item{\dots}{aditional arguments}
}

\details{generates a report of the queueing model received as parameter}

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

## Report the results
summary(o_MMInfKK)
}
\keyword{M/M/Inf/K/K}

