% File man/QueueingModel.i_MMInfKK.Rd
\name{QueueingModel.i_MMInfKK}
\alias{QueueingModel.i_MMInfKK}
\title{Builds a M/M/Inf/K/K queueing model}
\description{
  Builds a M/M/Inf/K/K queueing model
}
\usage{
  \method{QueueingModel}{i_MMInfKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMInfKK}
  \item{\dots}{aditional arguments}
}

\details{Build a M/M/Inf/K/K queueing model. It also checks the input params calling the \link{CheckInput.i_MMInfKK}}

\references{
  [Kleinrock1975] Leonard Kleinrock (1975).\cr
\emph{Queueing Systems Vol 1: Theory}.\cr
John Wiley & Sons.
}

\seealso{
  \code{\link{CheckInput.i_MMInfKK}}
}


\examples{
## create input parameters
i_MMInfKK <- NewInput.MMInfKK(lambda=0.25, mu=4, k=4)

## Build the model
QueueingModel(i_MMInfKK)

}
\keyword{M/M/Inf/K/K}
