% File man/QueueingModel.i_MMInf.Rd
\name{QueueingModel.i_MMInf}
\alias{QueueingModel.i_MMInf}
\title{Builds a M/M/Infinite queue model}
\description{
  Builds a M/M/Infinite queue model
}
\usage{
  \method{QueueingModel}{i_MMInf}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Build a M/M/Infinite model. It also checks the input params calling the \link{CheckInput.i_MMInf}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MMInf}}
}


\examples{
## create input parameters
i_mminf <- NewInput.MMInf(lambda=0.25, mu=4, n=0)

## Build the model
QueueingModel(i_mminf)

}
\keyword{M/M/Infinite}
