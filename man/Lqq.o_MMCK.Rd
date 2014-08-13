% File man/Lqq.o_MMCK.Rd
\name{Lqq.o_MMCK}
\alias{Lqq.o_MMCK}
\title{Returns the mean number of customers in queue when there is queue in the M/M/c/K queueing model}
\description{
  Returns the mean number of customers in queue when there is queue in the M/M/c/K queueing model
}
\usage{
  \method{Lqq}{o_MMCK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCK}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean number of customers in queue when there is queue in the M/M/c/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMCK}}.
}


\examples{
## See example 10.11 in reference [Sixto2004] for more details.
## create input parameters
i_mmck <- NewInput.MMCK(lambda=8, mu=4, c=5, k=12)

## Build the model
o_mmck <- QueueingModel(i_mmck)

## Returns the Lqq
Lqq(o_mmck)
}

\keyword{M/M/c/K}
