% File man/Inputs.o_MMC.Rd
\name{Inputs.o_MMC}
\alias{Inputs.o_MMC}
\title{Returns the input parameters of a M/M/c queueing model}
\description{
  Returns the inputs parameters of a already built M/M/c queueing model
}
\usage{
  \method{Inputs}{o_MMC}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMC}
  \item{\dots}{aditional arguments}
}

\details{Returns the input parameters of a M/M/c queueing model. The inputs parameters are created calling previously the \link{NewInput.MMC}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MMC}}.
}


\examples{
## See example 10.9 in reference [Sixto2004] for more details.
## create input parameters
i_mmc <- NewInput.MMC(lambda=5, mu=10, c=2, n=0, method=0)

## Build the model
o_mmc <- QueueingModel(i_mmc)

## Retunns the Inputs
Inputs(o_mmc)
}

\keyword{M/M/c}
