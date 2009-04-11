% File man/QueueingModel.i_MM1K.Rd
\name{QueueingModel.i_MM1K}
\alias{QueueingModel.i_MM1K}
\title{Builds a M/M/1/K queueing model}
\description{
  Builds a M/M/1/K queueing model
}
\usage{
  \method{QueueingModel}{i_MM1K}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MM1K}
  \item{\dots}{aditional arguments}
}

\details{Build a M/M/1/K queueing model. It also checks the input params calling the \link{CheckInput.i_MM1K}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MM1K}}.
}


\examples{
## See example 10.7 in reference [Sixto2004] for more details.
## create input parameters
i_mm1k <- NewInput.MM1K(lambda=5, mu=5.714, k=15)

## Build the model
QueueingModel(i_mm1k)

}
\keyword{M/M/1/K}

