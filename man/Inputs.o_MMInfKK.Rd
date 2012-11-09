% File man/Inputs.o_MMInfKK.Rd
\name{Inputs.o_MMInfKK}
\alias{Inputs.o_MMInfKK}
\title{Returns the input parameters of a M/M/Infinite/K/K queueing model}
\description{
  Returns the inputs parameters of a already built M/M/Infinite/K/K queueing model
}
\usage{
  \method{Inputs}{o_MMInfKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMInfKK}
  \item{\dots}{aditional arguments}
}

\details{Returns the input parameters of a M/M/Infinite/K/K queueing model. The inputs parameters are created calling previously the \link{NewInput.MMInfKK}}

\references{
  [Kleinrock1975] Leonard Kleinrock (1975).\cr
\emph{Queueing Systems Vol 1: Theory}.\cr
John Wiley & Sons.
}

\seealso{
  \code{\link{NewInput.MMInfKK}}.
}


\examples{
## create input parameters
i_MMInfKK <- NewInput.MMInfKK(lambda=0.25, mu=4, k=4)

## Build the model
o_MMInfKK <- QueueingModel(i_MMInfKK)

## Retunns the Inputs
Inputs(o_MMInfKK)
}

\keyword{M/M/Infinite/K/K}
