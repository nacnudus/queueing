% File man/Pn.o_MM1K.Rd
\name{Pn.o_MM1K}
\alias{Pn.o_MM1K}
\alias{Qn.o_MM1K}
\title{Returns the probabilities of a M/M/1/K queueing model}
\description{
  Pn returns the probabilities that a M/M/1/K queueing model has n customers.\cr
  Qn returns the probabilities that an arrival that enter the system see n customers in it.
}
\usage{
  \method{Pn}{o_MM1K}(x, \dots)
  \method{Qn}{o_MM1K}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1K}
  \item{\dots}{aditional arguments}
}

\details{Pn returns the probabilities that a M/M/1/K queueing model has n customers.\cr
Qn returns the probabilities that an arrival that enter the system see n customers.}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MM1K}}.
}


\examples{
## See example 10.7 in reference [Sixto2004] for more details.
## create input parameters
i_mm1k <- NewInput.MM1K(lambda=5, mu=5.714, k=15)

## Build the model
o_mm1k <- QueueingModel(i_mm1k)

## Returns the probabilities
Pn(o_mm1k)
Qn(o_mm1k)


}
\keyword{M/M/1/K}
