% File man/RO.Rd
\name{RO}
\alias{RO}
\title{Reports the server use of a queueing model}
\description{
  Reports the server use of a queueing model)
}
\usage{
 RO(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInfKK, o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Reports the server use of a queueing model (or network)}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{RO.o_MM1}}\cr
  \code{\link{RO.o_MMC}}\cr
  \code{\link{RO.o_MM1K}}\cr
  \code{\link{RO.o_MMCK}}\cr
  \code{\link{RO.o_MM1KK}}\cr
  \code{\link{RO.o_MMCKK}}\cr
  \code{\link{RO.o_MMCC}}\cr
  \code{\link{RO.o_MMCKM}}\cr
  \code{\link{RO.o_MMInfKK}}\cr
  \code{\link{RO.o_MMInf}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Report the use of the server
RO(o_mm1)
}

\keyword{Queueing Models}
