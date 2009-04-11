% File man/QueueingModel.i_MMCKM.Rd
\name{QueueingModel.i_MMCKM}
\alias{QueueingModel.i_MMCKM}
\title{Builds a M/M/c/K/m queueing model}
\description{
  Builds a M/M/c/K/m queueing model
}
\usage{
  \method{QueueingModel}{i_MMCKM}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMCKM}
  \item{\dots}{aditional arguments}
}

\details{Build a M/M/c/K/m queueing model. It also checks the input params calling the \link{CheckInput.i_MMCKM}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MMCKM}}
}


\examples{
## create input parameters
i_mmckm <- NewInput.MMCKM(lambda=0.25, mu=4, c=2, k=4, m=8)

## Build the model
QueueingModel(i_mmckm)

}
\keyword{M/M/c/K/m}
