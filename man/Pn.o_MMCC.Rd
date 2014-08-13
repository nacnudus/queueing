% File man/Pn.o_MMCC.Rd
\name{Pn.o_MMCC}
\alias{Pn.o_MMCC}
\alias{Qn.o_MMCC}
\title{Returns the probabilities of a M/M/c/c queueing model}
\description{
  Pn returns the probabilities that a M/M/c/c queueing model has n customers.\cr
  Qn returns the probabilities that an arrival that enter the system see n customers in it.
}
\usage{
  \method{Pn}{o_MMCC}(x, \dots)
  \method{Qn}{o_MMCC}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCC}
  \item{\dots}{aditional arguments}
}

\details{Pn returns the probabilities that a M/M/c/c queueing model has n customers.\cr
Qn returns the probabilities that an arrival that enter the system see n customers.}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMCC}}.
}

\examples{
## See example 10.12 in reference [Sixto2004] for more details.
## create input parameters
i_mmcc <- NewInput.MMCC(lambda=3, mu=0.25, c=15)

## Build the model
o_mmcc <- QueueingModel(i_mmcc)

## Returns the probabilities
Pn(o_mmcc)
Qn(o_mmcc)
}
\keyword{M/M/c/c}

