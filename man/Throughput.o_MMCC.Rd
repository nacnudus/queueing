% File man/Throughput.o_MMCC.Rd
\name{Throughput.o_MMCC}
\alias{Throughput.o_MMCC}
\title{Throughput of a M/M/c/c queueing model}
\description{
  Returns the throughput of a M/M/c/c queueing model
}
\usage{
  \method{Throughput}{o_MMCC}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCC}
  \item{\dots}{aditional arguments}
}

\details{Returns the throughput of a M/M/c/c queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MMCC}}, \code{\link{CheckInput.i_MMCC}}, \code{\link{QueueingModel.i_MMCC}}
}


\examples{
## See example 10.12 in reference [Sixto2004] for more details.
## create input parameters
i_mmcc <- NewInput.MMCC(lambda=3, mu=0.25, c=15)

## Build the model
o_mmcc <- QueueingModel(i_mmcc)

## Throughput
Throughput(o_mmcc)

}
\keyword{M/M/c/c}
