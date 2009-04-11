% File man/summary.o_MMCKK.Rd
\name{summary.o_MMCKK}
\alias{summary.o_MMCKK}
\title{Reports the results of a M/M/c/K/K queueing model}
\description{
  Reports the results of a M/M/c/K/K queueing model.
}
\usage{
  \method{summary}{o_MMCKK}(object, \dots)
}
\arguments{
  \item{object}{a object of class o_MMCKK}
  \item{\dots}{aditional arguments}
}

\details{generates a report of the queueing model received as parameter}

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

## Report the results
summary(o_mmckk)
}
\keyword{M/M/c/K/K}

