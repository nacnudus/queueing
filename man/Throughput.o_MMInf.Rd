% File man/Throughput.o_MMInf.Rd
\name{Throughput.o_MMInf}
\alias{Throughput.o_MMInf}
\title{Throughput of a M/M/Infinite queueing model}
\description{
  Returns the throughput of a M/M/Infinite queueing model
}
\usage{
  \method{Throughput}{o_MMInf}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Returns the throughput of a M/M/Infinite queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MMInf}}, \code{\link{CheckInput.i_MMInf}}, \code{\link{QueueingModel.i_MMInf}}
}


\examples{
## create input parameters
i_mminf <- NewInput.MMInf(lambda=0.25, mu=4, n=0)

## Build the model
o_mminf <- QueueingModel(i_mminf)

## Throughput
Throughput(o_mminf)

}
\keyword{M/M/Infinite}

