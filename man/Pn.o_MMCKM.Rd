% File man/Pn.o_MMCKM.Rd
\name{Pn.o_MMCKM}
\alias{Pn.o_MMCKM}
\alias{Qn.o_MMCKM}
\title{Returns the probabilities of a M/M/c/K/m queueing model}
\description{
  Pn returns the probabilities that a M/M/c/K/m queueing model has n customers.\cr
  Qn returns the probabilities that an arrival that enter the system see n customers in it.
}
\usage{
  \method{Pn}{o_MMCKM}(x, \dots)
  \method{Qn}{o_MMCKM}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCKM}
  \item{\dots}{aditional arguments}
}

\details{Pn returns the probabilities that a M/M/c/K/m queueing model has n customers.\cr
Qn returns the probabilities that an arrival that enter the system see n customers.}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMCKM}}.
}

\examples{
## create input parameters
i_mmckm <- NewInput.MMCKM(lambda=0.25, mu=4, c=2, k=4, m=8, method=0)

## Build the model
o_mmckm <- QueueingModel(i_mmckm)

## Returns the probabilities
Pn(o_mmckm)
Qn(o_mmckm)

}
\keyword{M/M/c/K/m}
