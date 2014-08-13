% File man/Pn.o_MM1.Rd
\name{Pn.o_MM1}
\alias{Pn.o_MM1}
\alias{Qn.o_MM1}
\title{Returns the probabilities of a M/M/1 queueing model}
\description{
  Pn returns the probabilities that a M/M/1 queueing model has n customers.\cr
  Qn returns the probabilities that an arrival that enter the system see n customers in it.
}
\usage{
  \method{Pn}{o_MM1}(x, \dots)
  \method{Qn}{o_MM1}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1}
  \item{\dots}{aditional arguments}
}

\details{Pn returns the probabilities that a M/M/1/K queueing model has n customers.\cr
Qn returns the probabilities that an arrival that enter the system see n customers. By the PASTA property, both probabilities has to be the same.}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MM1}}.
}


\examples{
## See example 10.3 in reference [Sixto2004] for more details.
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Returns the probabilities
Pn(o_mm1)
Qn(o_mm1)

}
\keyword{M/M/1}
