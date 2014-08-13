% File man/RO.o_MMInf.Rd
\name{RO.o_MMInf}
\alias{RO.o_MMInf}
\title{Reports the server use of a M/M/Infinite queueing model}
\description{
  Reports the server use of a M/M/Infinite queueing model
}
\usage{
  \method{RO}{o_MMInf}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Reports the server use of a M/M/Infinite queueing model. It should be noted that in this model, the RO parameter has a different meaning, its the traffic intensity and it coincides exactly with the average number of customers in the system (L)}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMInf}}\cr
  \code{\link{L.o_MMInf}}
}


\examples{
## create input parameters
i_mminf <- NewInput.MMInf(lambda=0.25, mu=4, n=0)

## Build the model
o_mminf <- QueueingModel(i_mminf)

## Report the use of the server
RO(o_mminf)
}
\keyword{M/M/Infinite}
