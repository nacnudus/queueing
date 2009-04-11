% File man/QueueingModel.i_MM1KK.Rd
\name{QueueingModel.i_MM1KK}
\alias{QueueingModel.i_MM1KK}
\title{Builds a M/M/1/K/K queueing model}
\description{
  Builds a M/M/1/K/K queueing model
}
\usage{
  \method{QueueingModel}{i_MM1KK}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{Build a M/M/1/K/K queueing model. It also checks the input params calling the \link{CheckInput.i_MM1KK}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MM1KK}}.
}


\examples{
## See example 10.13 in reference [Sixto2004] for more details.
## create input parameters
i_mm1kk <- NewInput.MM1KK(lambda=0.25, mu=4, k=2)

## Build the model
QueueingModel(i_mm1kk)

}
\keyword{M/M/1/K/K}
