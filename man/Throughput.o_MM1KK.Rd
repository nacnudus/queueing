% File man/Throughput.o_MM1KK.Rd
\name{Throughput.o_MM1KK}
\alias{Throughput.o_MM1KK}
\title{Throughput of a M/M/1/K/K queueing model}
\description{
  Returns the throughput of a M/M/1/K/K queueing model
}
\usage{
  \method{Throughput}{o_MM1KK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{Returns the throughput of a M/M/1/K/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MM1KK}}, \code{\link{CheckInput.i_MM1KK}}, \code{\link{QueueingModel.i_MM1KK}}
}


\examples{
## See example 10.13 in reference [Sixto2004] for more details.
## create input parameters
i_mm1kk <- NewInput.MM1KK(lambda=0.25, mu=4, k=2)

## Build the model
o_MM1KKk <- QueueingModel(i_mm1kk)

## Throughput
Throughput(o_MM1KKk)

}
\keyword{M/M/1/K/K}
