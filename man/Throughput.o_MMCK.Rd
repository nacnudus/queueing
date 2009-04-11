% File man/Throughput.o_MMCK.Rd
\name{Throughput.o_MMCK}
\alias{Throughput.o_MMCK}
\title{Throughput of a M/M/c/K queueing model}
\description{
  Returns the throughput of a M/M/c/K queueing model
}
\usage{
  \method{Throughput}{o_MMCK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCK}
  \item{\dots}{aditional arguments}
}

\details{Returns the throughput of a M/M/c/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MMCK}}, \code{\link{CheckInput.i_MMCK}}, \code{\link{QueueingModel.i_MMCK}}
}


\examples{
## See example 10.11 in reference [Sixto2004] for more details.
## create input parameters
i_mmck <- NewInput.MMCK(lambda=8, mu=4, c=5, k=12)

## Build the model
o_mmck <- QueueingModel(i_mmck)

## Throughput
Throughput(o_mmck)

}
\keyword{M/M/c/K}
