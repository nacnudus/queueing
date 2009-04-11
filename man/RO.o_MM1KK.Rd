% File man/RO.o_MM1KK.Rd
\name{RO.o_MM1KK}
\alias{RO.o_MM1KK}
\title{Reports the server use of a M/M/1/K/K queueing model}
\description{
  Reports the server use of a M/M/1/K/K queueing model
}
\usage{
  \method{RO}{o_MM1KK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{Reports the server use of a M/M/1/K/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MM1KK}}.
}


\examples{
## See example 10.13 in reference [Sixto2004] for more details.
## create input parameters
i_mm1kk <- NewInput.MM1KK(lambda=0.25, mu=4, k=2)

## Build the model
o_mm1kk <- QueueingModel(i_mm1kk)

## Report the use of the server
RO(o_mm1kk)
}
\keyword{M/M/1/K/K}

