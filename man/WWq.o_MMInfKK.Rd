% File man/WWq.o_MMInfKK.Rd
\name{WWq.o_MMInfKK}
\alias{WWq.o_MMInfKK}
\title{Returns the mean time spend in queue when there is queue in the M/M/Inf/K/K queueing model}
\description{
  Returns the mean time spend in queue when there is queue in the M/M/Inf/K/K queueing model
}
\usage{
  \method{WWq}{o_MMInfKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMInfKK}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean time spend in queue when there is queue in the M/M/Inf/K/K queueing model}

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

## Returns the WWq
WWq(o_MMInfKK)
}

\keyword{M/M/Inf/K/K}
