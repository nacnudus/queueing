% File man/Pn.o_MMCK.Rd
\name{Pn.o_MMCK}
\alias{Pn.o_MMCK}
\alias{Qn.o_MMCK}
\title{Returns the probabilities of a M/M/c/K queueing model}
\description{
  Pn returns the probabilities that a M/M/c/K queueing model has n customers.\cr
  Qn returns the probabilities that an arrival that enter the system see n customers in it.
}
\usage{
  \method{Pn}{o_MMCK}(x, \dots)
  \method{Qn}{o_MMCK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCK}
  \item{\dots}{aditional arguments}
}

\details{Pn returns the probabilities that a M/M/c/K queueing model has n customers.\cr
Qn returns the probabilities that an arrival that enter the system see n customers.}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMCK}}.
}

\examples{
## See example 10.11 in reference [Sixto2004] for more details.
## create input parameters
i_mmck <- NewInput.MMCK(lambda=8, mu=4, c=5, k=12)

## Build the model
o_mmck <- QueueingModel(i_mmck)

## Returns the probabilities
Pn(o_mmck)
Qn(o_mmck)
}
\keyword{M/M/c/K}

