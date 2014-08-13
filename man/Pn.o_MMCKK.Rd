% File man/Pn.o_MMCKK.Rd
\name{Pn.o_MMCKK}
\alias{Pn.o_MMCKK}
\alias{Qn.o_MMCKK}
\title{Returns the probabilities of a M/M/c/K/K queueing model}
\description{
  Pn returns the probabilities that a M/M/c/K/K queueing model has n customers.\cr
  Qn returns the probabilities that an arrival that enter the system see n customers in it.
}
\usage{
  \method{Pn}{o_MMCKK}(x, \dots)
  \method{Qn}{o_MMCKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCKK}
  \item{\dots}{aditional arguments}
}

\details{Pn returns the probabilities that a M/M/c/K/K queueing model has n customers.\cr
Qn returns the probabilities that an arrival that enter the system see n customers.}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMCKK}}.
}


\examples{
## create input parameters
i_mmckk <- NewInput.MMCKK(lambda=8, mu=2, c=5, k=12, method=0)

## Build the model
o_mmckk <- QueueingModel(i_mmckk)

## Returns the parameters
Pn(o_mmckk)
Qn(o_mmckk)
}
\keyword{M/M/c/K/K}
