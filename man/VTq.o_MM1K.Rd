% File man/VTq.o_MM1K.Rd
\name{VTq.o_MM1K}
\alias{VTq.o_MM1K}
\title{Returns the variance of the time spend in queue in the M/M/1/K queueing model}
\description{
  Returns the variance of the time spend in queue in the M/M/1/K queueing model
}
\usage{
  \method{VTq}{o_MM1K}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1K}
  \item{\dots}{aditional arguments}
}

\details{Returns the variance of the time spend in queue in the M/M/1/K queueing model}

\references{
  [Sztrik2012] Dr. Janos Sztrik (2012).\cr
\emph{Basic Queueing Theory}.\cr
University of Debrecen, Faculty of Informatics.
}

\seealso{
  \code{\link{QueueingModel.i_MM1K}}.
}


\examples{

## create input parameters
i_mm1k <- NewInput.MM1K(lambda=5, mu=5.714, k=15)

## Build the model
o_mm1k <- QueueingModel(i_mm1k)

## Returns the variance
VTq(o_mm1k)
}

\keyword{M/M/1/K}
