% File man/RO.o_MMCC.Rd
\name{RO.o_MMCC}
\alias{RO.o_MMCC}
\title{Reports the server use of a M/M/c/c queueing model}
\description{
  Reports the server use of a M/M/c/c queueing model
}
\usage{
  \method{RO}{o_MMCC}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCC}
  \item{\dots}{aditional arguments}
}

\details{Reports the server use of a M/M/c/c queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMCC}}.
}


\examples{
## See example 10.12 in reference [Sixto2004] for more details.
## create input parameters
i_mmcc <- NewInput.MMCC(lambda=3, mu=0.25, c=15)

## Build the model
o_mmcc <- QueueingModel(i_mmcc)

## Report the use of the server
RO(o_mmcc)
}
\keyword{M/M/c/c}

