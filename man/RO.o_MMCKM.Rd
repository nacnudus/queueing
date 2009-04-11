% File man/RO.o_MMCKM.Rd
\name{RO.o_MMCKM}
\alias{RO.o_MMCKM}
\title{Reports the server use of a M/M/c/K/m queueing model}
\description{
  Reports the server use of a M/M/c/K/m queueing model
}
\usage{
  \method{RO}{o_MMCKM}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCKM}
  \item{\dots}{aditional arguments}
}

\details{Reports the server use of a M/M/c/K/m queueing model}

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

## Report the use of the server
RO(o_mmckm)
}
\keyword{M/M/c/K/m}
