% File man/VN.Rd
\name{VN}
\alias{VN}
\title{Returns the variance of the number of customers in a queueing model (or network)}
\description{
  Returns the variance of the number of customers in a queueing model (or network)
}
\usage{
  VN(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInfKK, o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Returns the variance of the number of customers in a queueing model (or network)}

\references{
  [Sztrik2012] Dr. Janos Sztrik (2012).\cr
\emph{Basic Queueing Theory}.\cr
University of Debrecen, Faculty of Informatics.
}

\seealso{
  \code{\link{VN.o_MM1}}\cr
  \code{\link{VN.o_MMC}}\cr
  \code{\link{VN.o_MMCC}}\cr
  \code{\link{VN.o_MMInf}}\cr
  \code{\link{VN.o_MMInfKK}}\cr
  \code{\link{VN.o_MM1K}}\cr
  \code{\link{VN.o_MMCK}}\cr
  \code{\link{VN.o_MM1KK}}\cr
  \code{\link{VN.o_MMCKK}}\cr
  \code{\link{VN.o_MMCKM}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Returns the variance
VN(o_mm1)
}

\keyword{Queueing Models}
