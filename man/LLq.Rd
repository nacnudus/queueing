% File man/LLq.Rd
\name{LLq}
\alias{LLq}
\title{Returns the mean number of clients in queue when there is queue in a queueing model}
\description{
  Returns the mean number of clients in queue when there is queue in a queueing model
}
\usage{
  LLq(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1, o_MMInfKK, o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean number of clients in queue when there is queue in a queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{LLq.o_MM1}}\cr
  \code{\link{LLq.o_MMInfKK}}\cr
  \code{\link{LLq.o_MMInf}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Returns the LLq
LLq(o_mm1)
}

\keyword{M/M/1}
