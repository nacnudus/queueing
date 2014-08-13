% File man/VTq.Rd
\name{VTq}
\alias{VTq}
\title{Returns the variance of the time spend in queue in a queueing model}
\description{
  Returns the variance of the time spend in queue in a queueing model
}
\usage{
  VTq(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInfKK, o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Returns the variance of the time spend in queue in a queueing model}

\references{
  [Sztrik2012] Dr. Janos Sztrik (2012).\cr
\emph{Basic Queueing Theory}.\cr
University of Debrecen, Faculty of Informatics.
}

\seealso{
  \code{\link{VTq.o_MM1}}\cr
  \code{\link{VTq.o_MMC}}\cr
  \code{\link{VTq.o_MMCC}}\cr
  \code{\link{VTq.o_MMInf}}\cr
  \code{\link{VTq.o_MMInfKK}}\cr
  \code{\link{VTq.o_MM1K}}\cr
  \code{\link{VTq.o_MMCK}}\cr
  \code{\link{VTq.o_MM1KK}}\cr
  \code{\link{VTq.o_MMCKK}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Returns the variance of the time spend in queue
VTq(o_mm1)
}

\keyword{M/M/1}
