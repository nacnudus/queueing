% File man/RO.o_MMC.Rd
\name{RO.o_MMC}
\alias{RO.o_MMC}
\title{Reports the server use of a M/M/c queueing model}
\description{
  Reports the server use of a M/M/c queueing model
}
\usage{
  \method{RO}{o_MMC}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMC}
  \item{\dots}{aditional arguments}
}

\details{Reports the server use of a M/M/c queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMC}}.
}


\examples{
## See example 10.9 in reference [Sixto2004] for more details.
## create input parameters
i_mmc <- NewInput.MMC(lambda=5, mu=10, c=2, n=0)

## Build the model
o_mmc <- QueueingModel(i_mmc)

## Report the use of the server
RO(o_mmc)
}
\keyword{M/M/c}
