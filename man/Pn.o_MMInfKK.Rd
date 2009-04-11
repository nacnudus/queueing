% File man/Pn.o_MMInfKK.Rd
\name{Pn.o_MMInfKK}
\alias{Pn.o_MMInfKK}
\title{Returns the probabilities of a M/M/Inf/K/K queueing model}
\description{
  Returns the probabilities of a M/M/Inf/K/K queueing model
}
\usage{
  \method{Pn}{o_MMInfKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMInfKK}
  \item{\dots}{aditional arguments}
}

\details{Returns the probabilities of a M/M/Inf/K/K queueing model}

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

## Returns the probabilities
Pn(o_MMInfKK)

}
\keyword{M/M/Inf/K/K}
