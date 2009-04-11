% File man/Throughput.o_MM1.Rd
\name{Throughput.o_MM1}
\alias{Throughput.o_MM1}
\title{Throughput of a M/M/1 queueing model}
\description{
  Returns the throughput of a M/M/1 queueing model
}
\usage{
  \method{Throughput}{o_MM1}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1}
  \item{\dots}{aditional arguments}
}

\details{Returns the throughput of a M/M/1 queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MM1}}, \code{\link{CheckInput.i_MM1}}, \code{\link{QueueingModel.i_MM1}}
}


\examples{
## See example 10.3 in reference [Sixto2004] for more details.
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Throughput
Throughput(o_mm1)

}
\keyword{M/M/1}
