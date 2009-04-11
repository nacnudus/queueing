% File man/summary.o_MMCKM.Rd
\name{summary.o_MMCKM}
\alias{summary.o_MMCKM}
\title{Reports the results of a M/M/c/K/m queueing model}
\description{
  Reports the results of a M/M/c/K/m queueing model.
}
\usage{
  \method{summary}{o_MMCKM}(object, \dots)
}
\arguments{
  \item{object}{a object of class o_MMCKM}
  \item{\dots}{aditional arguments}
}

\details{generates a report of the queueing model received as parameter}

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
i_mmckm <- NewInput.MMCKM(lambda=0.25, mu=4, c=2, k=4, m=8)

## Build the model
o_mmckm <- QueueingModel(i_mmckm)

## Report the results
summary(o_mmckm)
}
\keyword{M/M/c/K/m}

