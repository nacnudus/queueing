% File man/Pn.o_MMInf.Rd
\name{Pn.o_MMInf}
\alias{Pn.o_MMInf}
\alias{Qn.o_MMInf}
\title{Returns the probabilities of a M/M/Infinite queueing model}
\description{
  Pn returns the probabilities that a M/M/Infinite queueing model has n clients.\cr
  Qn returns the probabilities that an arrival that enter the system see n clients in it.
}
\usage{
  \method{Pn}{o_MMInf}(x, \dots)
  \method{Qn}{o_MMInf}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Pn returns the probabilities that a M/M/Infinite queueing model has n clients.\cr
Qn returns the probabilities that an arrival that enter the system see n clients. By the PASTA property, both probabilities has to be the same.}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMInf}}.
}

\examples{
## create input parameters
i_mminf <- NewInput.MMInf(lambda=0.25, mu=4, n=0)

## Build the model
o_mminf <- QueueingModel(i_mminf)

## Returns the probabilities
Pn(o_mminf)
Qn(o_mminf)

}
\keyword{M/M/Infinite}
