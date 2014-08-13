% File man/VT.Rd
\name{VT}
\alias{VT}
\title{Returns the variance of the time spend in a queueing model (or network)}
\description{
  Returns the variance of the time spend in a queueing model (or network)
}
\usage{
  VT(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInfKK, o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Returns the variance of the time spend in a queueing model (or network)}

\references{
  [Sztrik2012] Dr. Janos Sztrik (2012).\cr
\emph{Basic Queueing Theory}.\cr
University of Debrecen, Faculty of Informatics.
}

\seealso{
  \code{\link{VT.o_MM1}}\cr
  \code{\link{VT.o_MMC}}\cr
  \code{\link{VT.o_MMCC}}\cr
  \code{\link{VT.o_MMInf}}\cr
  \code{\link{VT.o_MMInfKK}}\cr
  \code{\link{VT.o_MM1K}}\cr
  \code{\link{VT.o_MM1KK}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Returns the variance of the time spend in the system
VT(o_mm1)
}

\keyword{Queueing Models}
