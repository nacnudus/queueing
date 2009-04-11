% File man/Inputs.Rd
\name{Inputs}
\alias{Inputs}
\title{Returns the input parameters of a queueing model (or network)}
\description{
  Returns the inputs parameters of a already built queueing model (or network)
}
\usage{
  Inputs(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInfKK, o_MMInf, o_OJN, o_MCON}
  \item{\dots}{aditional arguments}
}

\details{Returns the input parameters of a queueing model (or network)}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{Inputs.o_MM1}}\cr
  \code{\link{Inputs.o_MMC}}\cr
  \code{\link{Inputs.o_MM1K}}\cr
  \code{\link{Inputs.o_MMCK}}\cr
  \code{\link{Inputs.o_MM1KK}}\cr
  \code{\link{Inputs.o_MMCKK}}\cr
  \code{\link{Inputs.o_MMCC}}\cr
  \code{\link{Inputs.o_MMCKM}}\cr
  \code{\link{Inputs.o_MMInfKK}}\cr
  \code{\link{Inputs.o_MMInf}}\cr
  \code{\link{Inputs.o_OJN}}\cr
  \code{\link{Inputs.o_MCON}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Retunns the Inputs
Inputs(o_mm1)
}

\keyword{Queueing Models}
