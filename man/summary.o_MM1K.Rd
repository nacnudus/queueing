% File man/summary.o_MM1K.Rd
\name{summary.o_MM1K}
\alias{summary.o_MM1K}
\title{Reports the results of a M/M/1/K queueing model}
\description{
  Reports the results of a M/M/1/K queueing model.
}
\usage{
  \method{summary}{o_MM1K}(object, \dots)
}
\arguments{
  \item{object}{a object of class o_MM1K}
  \item{\dots}{aditional arguments}
}

\details{generates a report of the queueing model received as parameter}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MM1K}}.
}


\examples{
## See example 10.7 in reference [Sixto2004] for more details.
## create input parameters
i_mm1k <- NewInput.MM1K(lambda=5, mu=5.714, k=15)

## Build the model
o_mm1k <- QueueingModel(i_mm1k)

## Report the results
summary(o_mm1k)
}
\keyword{M/M/1/K}
