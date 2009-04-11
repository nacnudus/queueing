% File man/Inputs.o_MMCC.Rd
\name{Inputs.o_MMCC}
\alias{Inputs.o_MMCC}
\title{Returns the input parameters of a M/M/c/c queueing model}
\description{
  Returns the inputs parameters of a already built M/M/c/c queueing model
}
\usage{
  \method{Inputs}{o_MMCC}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCC}
  \item{\dots}{aditional arguments}
}

\details{Returns the input parameters of a M/M/c/c queueing model. The inputs parameters are created calling previously the \link{NewInput.MMCC}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MMCC}}.
}


\examples{
## See example 10.12 in reference [Sixto2004] for more details.
## create input parameters
i_mmcc <- NewInput.MMCC(lambda=3, mu=0.25, c=15)

## Build the model
o_mmcc <- QueueingModel(i_mmcc)

## Retunns the Inputs
Inputs(o_mmcc)
}

\keyword{M/M/c/c}
