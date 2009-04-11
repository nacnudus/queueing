% File man/W.o_MMCKM.Rd
\name{W.o_MMCKM}
\alias{W.o_MMCKM}
\title{Returns the mean time spend in the M/M/c/K/m queueing model}
\description{
  Returns the mean time spend in the M/M/c/K/m queueing model
}
\usage{
  \method{W}{o_MMCKM}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCKM}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean time spend in the M/M/c/K/m queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMCKM}}.
}


\examples{
## create input parameters
i_mmckm <- NewInput.MMCKM(lambda=0.25, mu=4, c=2, k=4, m=8, method=0)

## Build the model
o_mmckm <- QueueingModel(i_mmckm)

## Returns the W
W(o_mmckm)
}

\keyword{M/M/c/K/m}
