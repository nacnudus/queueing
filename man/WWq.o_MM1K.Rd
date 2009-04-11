% File man/WWq.o_MM1K.Rd
\name{WWq.o_MM1K}
\alias{WWq.o_MM1K}
\title{Returns the mean time spend in queue when there is queue in the M/M/1/K queueing model}
\description{
  Returns the mean time spend in queue when there is queue in the M/M/1/K queueing model
}
\usage{
  \method{WWq}{o_MM1K}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1K}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean time spend in queue when there is queue in the M/M/1/K queueing model}

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

## Returns the WWq
WWq(o_mm1k)
}

\keyword{M/M/1/K}
