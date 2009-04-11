% File man/QueueingModel.i_MMCKK.Rd
\name{QueueingModel.i_MMCKK}
\alias{QueueingModel.i_MMCKK}
\title{Builds a M/M/c/K/K queueing model}
\description{
  Builds a M/M/c/K/K queueing model
}
\usage{
  \method{QueueingModel}{i_MMCKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMCKK}
  \item{\dots}{aditional arguments}
}

\details{Build a M/M/c/K/K queueing model. It also checks the input params calling the \link{CheckInput.i_MMCKK}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MMCKK}}.
}


\examples{
## create input parameters
i_mmckk <- NewInput.MMCKK(lambda=8, mu=2, c=5, k=12)

## Build the model
QueueingModel(i_mmckk)

}
\keyword{M/M/c/K/K}
