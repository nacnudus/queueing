% File man/QueueingModel.i_MM1.Rd
\name{QueueingModel.i_MM1}
\alias{QueueingModel.i_MM1}
\title{Builds a M/M/1 queueing model}
\description{
  Builds a M/M/1 queueing model
}
\usage{
  \method{QueueingModel}{i_MM1}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MM1}
  \item{\dots}{aditional arguments}
}

\details{Build a M/M/1 queueing model. It also checks the input params calling the \link{CheckInput.i_MM1}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MM1}}
}


\examples{
## See example 10.3 in reference [Sixto2004] for more details.
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
QueueingModel(i_mm1)

}
\keyword{M/M/1}
