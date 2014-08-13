% File man/VTq.o_MM1KK.Rd
\name{VTq.o_MM1KK}
\alias{VTq.o_MM1KK}
\title{Returns the variance of the time spend in queue in the M/M/1/K/K queueing model}
\description{
  Returns the variance of the time spend in queue in the M/M/1/K/K queueing model
}
\usage{
  \method{VTq}{o_MM1KK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{Returns the variance of the time spend in queue in the M/M/1/K/K queueing model}

\references{
  [Sztrik2012] Dr. Janos Sztrik (2012).\cr
\emph{Basic Queueing Theory}.\cr
University of Debrecen, Faculty of Informatics.
}

\seealso{
  \code{\link{QueueingModel.i_MM1KK}}.
}


\examples{
## See example 10.13 in reference [Sixto2004] for more details.
## create input parameters
i_mm1kk <- NewInput.MM1KK(lambda=0.25, mu=4, k=2, method=3)

## Build the model
o_mm1kk <- QueueingModel(i_mm1kk)

## Returns the VTq
VTq(o_mm1kk)
}

\keyword{M/M/1/K/K}
