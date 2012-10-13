% File man/W.Rd
\name{W}
\alias{W}
\title{Returns the mean time spend in a queueing model (or network)}
\description{
  Returns the mean time spend in a queueing model (or network)
}
\usage{
  W(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInfKK, o_MMInf, o_OJN, o_MCON, o_MCCN}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean time spend in a queueing model (or network)}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{W.o_MM1}}\cr
  \code{\link{W.o_MMC}}\cr
  \code{\link{W.o_MM1K}}\cr
  \code{\link{W.o_MMCK}}\cr
  \code{\link{W.o_MM1KK}}\cr
  \code{\link{W.o_MMCKK}}\cr
  \code{\link{W.o_MMCC}}\cr
  \code{\link{W.o_MMCKM}}\cr
  \code{\link{W.o_MMInfKK}}\cr
  \code{\link{W.o_MMInf}}\cr
  \code{\link{W.o_OJN}}\cr
  \code{\link{W.o_MCON}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Returns the W
W(o_mm1)
}

\keyword{Queueing Models}
