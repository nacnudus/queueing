% File man/QueueingModel.i_MMC.Rd
\name{QueueingModel.i_MMC}
\alias{QueueingModel.i_MMC}
\title{Builds a M/M/c queueing model}
\description{
  Builds a M/M/c queueing model
}
\usage{
  \method{QueueingModel}{i_MMC}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMC}
  \item{\dots}{aditional arguments}
}

\details{Build a M/M/c/ queueing model. It also checks the input params calling the \link{CheckInput.i_MMC}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MMC}}
}


\examples{
## See example 10.9 in reference [Sixto2004] for more details.
## create input parameters
i_mmc <- NewInput.MMC(lambda=5, mu=10, c=2, n=0, method=0)

## Build the model
QueueingModel(i_mmc)

}
\keyword{M/M/c}

