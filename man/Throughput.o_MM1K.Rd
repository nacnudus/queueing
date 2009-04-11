% File man/Throughput.o_MM1K.Rd
\name{Throughput.o_MM1K}
\alias{Throughput.o_MM1K}
\title{Throughput of a M/M/1/K queueing model}
\description{
  Returns the throughput of a M/M/1/K queueing model
}
\usage{
  \method{Throughput}{o_MM1K}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1K}
  \item{\dots}{aditional arguments}
}

\details{Returns the throughput of a M/M/1/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MM1K}}, \code{\link{CheckInput.i_MM1K}}, \code{\link{QueueingModel.i_MM1K}}
}


\examples{
## See example 10.7 in reference [Sixto2004] for more details.
## create input parameters
i_mm1k <- NewInput.MM1K(lambda=5, mu=5.714, k=15)

## Build the model
o_mmck <- QueueingModel(i_mm1k)

## Throughput
Throughput(o_mmck)

}
\keyword{M/M/1/K}

