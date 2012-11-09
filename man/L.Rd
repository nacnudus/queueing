% File man/L.Rd
\name{L}
\alias{L}
\title{Returns the mean number of clients in a queueing model (or network)}
\description{
  Returns the mean number of clients in a queueing model (or network)
}
\usage{
  L(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInfKK, o_MMInf, o_OJN, o_MCON, o_MCCN, o_MCMN}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean number of clients in a queueing model (or network)}

\references{
 [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{L.o_MM1}}\cr
  \code{\link{L.o_MMC}}\cr
  \code{\link{L.o_MM1K}}\cr
  \code{\link{L.o_MMCK}}\cr
  \code{\link{L.o_MM1KK}}\cr
  \code{\link{L.o_MMCKK}}\cr
  \code{\link{L.o_MMCC}}\cr
  \code{\link{L.o_MMCKM}}\cr
  \code{\link{L.o_MMInfKK}}\cr
  \code{\link{L.o_MMInf}}\cr
  \code{\link{L.o_OJN}}\cr
  \code{\link{L.o_CJN}}\cr
  \code{\link{L.o_MCON}}\cr
  \code{\link{L.o_MCCN}}\cr
  \code{\link{L.o_MCMN}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Retunns the L
L(o_mm1)
}

\keyword{Queueing Models}
