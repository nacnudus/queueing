% File man/L.o_MM1KK.Rd
\name{L.o_MM1KK}
\alias{L.o_MM1KK}
\title{Returns the mean number of customers in the M/M/1/K/K queueing model}
\description{
  Returns the mean number of customers in the M/M/1/K/K queueing model
}
\usage{
  \method{L}{o_MM1KK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean number of customers in the M/M/1/K/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MM1K}}.
}


\examples{
## See example 10.13 in reference [Sixto2004] for more details.
## create input parameters
i_mm1kk <- NewInput.MM1KK(lambda=0.25, mu=4, k=2, method=3)

## Build the model
o_mm1kk <- QueueingModel(i_mm1kk)

## Returns the L
L(o_mm1kk)
}

\keyword{M/M/1/K/K}
