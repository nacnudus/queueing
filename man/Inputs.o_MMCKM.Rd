% File man/Inputs.o_MMCKM.Rd
\name{Inputs.o_MMCKM}
\alias{Inputs.o_MMCKM}
\title{Returns the input parameters of a M/M/c/K/m queueing model}
\description{
  Returns the inputs parameters of a already built M/M/c/K/m queueing model
}
\usage{
  \method{Inputs}{o_MMCKM}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCKM}
  \item{\dots}{aditional arguments}
}

\details{Returns the input parameters of a M/M/c/K/m queueing model. The inputs parameters are created calling previously the \link{NewInput.MMCKM}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MMCKM}}.
}


\examples{
## create input parameters
i_mmckm <- NewInput.MMCKM(lambda=0.25, mu=4, c=2, k=4, m=8)

## Build the model
o_mmckm <- QueueingModel(i_mmckm)

## Retunns the Inputs
Inputs(o_mmckm)
}

\keyword{M/M/c/K/m}
