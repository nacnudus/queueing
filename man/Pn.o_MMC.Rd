% File man/Pn.o_MMC.Rd
\name{Pn.o_MMC}
\alias{Pn.o_MMC}
\alias{Qn.o_MMC}
\title{Returns the probabilities of a M/M/c queueing model}
\description{
  Pn returns the probabilities that a M/M/c queueing model has n clients.\cr
  Qn returns the probabilities that an arrival that enter the system see n clients in it.
}
\usage{
  \method{Pn}{o_MMC}(x, \dots)
  \method{Qn}{o_MMC}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMC}
  \item{\dots}{aditional arguments}
}

\details{Pn returns the probabilities that a M/M/c queueing model has n clients.\cr
Qn returns the probabilities that an arrival that enter the system see n clients. By the PASTA property, both probabilities has to be the same.}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMC}}.
}


\examples{
## See example 10.9 in reference [Sixto2004] for more details.
## create input parameters
i_mmc <- NewInput.MMC(lambda=5, mu=10, c=2, n=0, method=0)

## Build the model
o_mmc <- QueueingModel(i_mmc)

## Returns the probabilities
Pn(o_mmc)
Qn(o_mmc)
}
\keyword{M/M/c}

