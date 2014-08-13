% File man/SP.o_MM1KK.Rd
\name{SP.o_MM1KK}
\alias{SP.o_MM1KK}
\title{Returns the saturation point of a M/M/1/K/K queueing model}
\description{
  Returns the saturation point, or the maximum number of customers that the M/M/1/K/K queueing model can support with no interference or syncronization between themselves
}
\usage{
  \method{SP}{o_MM1KK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{The value returned is the optimal number of customers of a M/M/1/K/K queueing model. It coincides with the inverse of the serialization parameter of Amdahl's Law. That is, the value which converges the speedup func(k) = k/(1 + ser * (k-1)). It makes sense, because the saturation point is the maximun value in which no syncronization happens.}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MM1KK}}
}


\examples{
## See example 10.13 in reference [Sixto2004] for more details.
## create input parameters
i_mm1kk <- NewInput.MM1KK(lambda=0.25, mu=4, k=4, method=3)

## Build the model
o_mm1kk <- QueueingModel(i_mm1kk)

## Returns the saturation point
SP(o_mm1kk)

}
\keyword{M/M/1/K/K}
