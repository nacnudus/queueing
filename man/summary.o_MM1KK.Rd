% File man/summary.o_MM1KK.Rd
\name{summary.o_MM1KK}
\alias{summary.o_MM1KK}
\title{Reports the results of a M/M/1/K/K queueing model}
\description{
  Reports the results of a M/M/1/K/K queueing model.
}
\usage{
  \method{summary}{o_MM1KK}(object, \dots)
}
\arguments{
  \item{object}{a object of class o_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{generates a report of the queueing model received as parameter}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MM1KK}}.
}


\examples{
## See example 10.13 in reference [Sixto2004] for more details.
## create input parameters
i_mm1kk <- NewInput.MM1KK(lambda=0.25, mu=4, k=2, method=0)

## Build the model
o_mm1kk <- QueueingModel(i_mm1kk)

## Report the results
summary(o_mm1kk)
}
\keyword{M/M/1/K/K}
