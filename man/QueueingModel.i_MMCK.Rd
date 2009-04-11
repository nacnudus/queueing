% File man/QueueingModel.i_MMCK.Rd
\name{QueueingModel.i_MMCK}
\alias{QueueingModel.i_MMCK}
\title{Builds a M/M/c/K queueing model}
\description{
  Builds a M/M/c/K queueing model
}
\usage{
  \method{QueueingModel}{i_MMCK}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMCK}
  \item{\dots}{aditional arguments}
}

\details{Build a M/M/c/K queueing model. It also checks the input params calling the \link{CheckInput.i_MMCK}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MMCK}}.
}


\examples{
## See example 10.11 in reference [Sixto2004] for more details.
## create input parameters
i_mmck <- NewInput.MMCK(lambda=8, mu=4, c=5, k=12)

## Build the model
QueueingModel(i_mmck)

}
\keyword{M/M/c/K}

