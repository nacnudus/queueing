% File man/Wq.o_MMInf.Rd
\name{Wq.o_MMInf}
\alias{Wq.o_MMInf}
\title{Returns the mean time spend in queue in the M/M/Infinite queueing model}
\description{
  Returns the mean time spend in queue in the M/M/Infinite queueing model
}
\usage{
  \method{Wq}{o_MMInf}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean time spend in queue in the M/M/Infinite queueing model}

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

## Returns the Wq
Wq(o_mminf)
}

\keyword{M/M/Infinite}
