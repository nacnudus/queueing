% File man/VTq.o_MMCKK.Rd
\name{VTq.o_MMCKK}
\alias{VTq.o_MMCKK}
\title{Returns the variance of the time spend in queue in the M/M/c/K/K queueing model}
\description{
  Returns the variance of the time spend in queue in the M/M/c/K/K queueing model
}
\usage{
  \method{VTq}{o_MMCKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCKK}
  \item{\dots}{aditional arguments}
}

\details{Returns the variance of the time spend in queue in the M/M/c/K/K queueing model}

\references{
  [Sztrik2012] Dr. Janos Sztrik (2012).\cr
\emph{Basic Queueing Theory}.\cr
University of Debrecen, Faculty of Informatics.
}

\seealso{
  \code{\link{QueueingModel.i_MMCKK}}.
}


\examples{
## create input parameters
i_mmckk <- NewInput.MMCKK(lambda=8, mu=2, c=5, k=12, method=0)

## Build the model
o_mmckk <- QueueingModel(i_mmckk)

## Returns the variance
VTq(o_mmckk)
}

\keyword{M/M/c/K/K}
