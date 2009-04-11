% File man/Inputs.o_MMCKK.Rd
\name{Inputs.o_MMCKK}
\alias{Inputs.o_MMCKK}
\title{Returns the input parameters of a M/M/c/K/K queueing model}
\description{
  Returns the inputs parameters of a already built M/M/c/K/K queueing model
}
\usage{
  \method{Inputs}{o_MMCKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCKK}
  \item{\dots}{aditional arguments}
}

\details{Returns the input parameters of a M/M/c/K/K queueing model. The inputs parameters are created calling previously the \link{NewInput.MMCKK}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MMCKK}}.
}


\examples{
## create input parameters
i_mmckk <- NewInput.MMCKK(lambda=8, mu=2, c=5, k=12)

## Build the model
o_mmckk <- QueueingModel(i_mmckk)

## Retunns the Inputs
Inputs(o_mmckk)
}

\keyword{M/M/c/K/K}
