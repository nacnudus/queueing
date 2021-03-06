% File man/Inputs.o_MMCK.Rd
\name{Inputs.o_MMCK}
\alias{Inputs.o_MMCK}
\title{Returns the input parameters of a M/M/c/K queueing model}
\description{
  Returns the inputs parameters of a already built M/M/c/K queueing model
}
\usage{
  \method{Inputs}{o_MMCK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCK}
  \item{\dots}{aditional arguments}
}

\details{Returns the input parameters of a M/M/c/K queueing model. The inputs parameters are created calling previously the \link{NewInput.MMCK}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MMCK}}.
}


\examples{
## See example 10.11 in reference [Sixto2004] for more details.
## create input parameters
i_mmck <- NewInput.MMCK(lambda=8, mu=4, c=5, k=12)

## Build the model
o_mmck <- QueueingModel(i_mmck)

## Retunns the Inputs
Inputs(o_mmck)
}

\keyword{M/M/c/K}
