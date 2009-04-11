% File man/summary.o_MMInf.Rd
\name{summary.o_MMInf}
\alias{summary.o_MMInf}
\title{Reports the results of a M/M/Infinite queueing model}
\description{
  Reports the results of a M/M/Infinite queueing model.
}
\usage{
  \method{summary}{o_MMInf}(object, \dots)
}
\arguments{
  \item{object}{a object of class o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{generates a report of the queueing model received as parameter}

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

## Report the results
summary(o_mminf)
}
\keyword{M/M/Infinite}
