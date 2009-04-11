% File man/QueueingModel.i_MMCC.Rd
\name{QueueingModel.i_MMCC}
\alias{QueueingModel.i_MMCC}
\title{Builds a M/M/c/c queueing model}
\description{
  Builds a M/M/c/c queueing model
}
\usage{
  \method{QueueingModel}{i_MMCC}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMCC}
  \item{\dots}{aditional arguments}
}

\details{Build a M/M/c/c queueing model. It also checks the input params calling the \link{CheckInput.i_MMCC}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MMCC}}.
}


\examples{
## See example 10.12 in reference [Sixto2004] for more details.
## create input parameters
i_mmcc <- NewInput.MMCC(lambda=3, mu=0.25, c=15)

## Build the model
QueueingModel(i_mmcc)

}
\keyword{M/M/c/c}

