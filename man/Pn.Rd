% File man/Pn.Rd
\name{Pn}
\alias{Pn}
\alias{Qn}
\title{Returns the probabilities of a queueing model (or network)}
\description{
  Pn returns the probabilities that a queueing model (or network) has n customers.\cr
  Qn returns the probabilities that an arrival that enter the system see n customers in it
}
\usage{
  Pn(x, \dots)
  Qn(x, \dots)
}
\arguments{
  \item{x}{For Pn, an object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInfKK, o_MMInf, o_OJN. For Qn, an object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInfKK, o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Pn returns the system probabilities of a queueing model (or network). Qn returns the probability that an effective arrival see n customers in the system}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{Pn.o_MM1}}\cr
  \code{\link{Qn.o_MM1}}\cr
  \code{\link{Pn.o_MMC}}\cr
  \code{\link{Qn.o_MMC}}\cr
  \code{\link{Pn.o_MM1K}}\cr
  \code{\link{Qn.o_MM1K}}\cr
  \code{\link{Pn.o_MMCK}}\cr
  \code{\link{Qn.o_MMCK}}\cr
  \code{\link{Pn.o_MM1KK}}\cr
  \code{\link{Qn.o_MM1KK}}\cr
  \code{\link{Pn.o_MMCKK}}\cr
  \code{\link{Qn.o_MMCKK}}\cr
  \code{\link{Pn.o_MMCC}}\cr
  \code{\link{Qn.o_MMCC}}\cr
  \code{\link{Pn.o_MMCKM}}\cr
  \code{\link{Qn.o_MMCKM}}\cr
  \code{\link{Pn.o_MMInfKK}}\cr
  \code{\link{Qn.o_MMInfKK}}\cr
  \code{\link{Pn.o_MMInf}}\cr
  \code{\link{Qn.o_MMInf}}\cr
  \code{\link{Pn.o_OJN}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Returns the probabilities
Pn(o_mm1)

}
\keyword{Queueing Models}
