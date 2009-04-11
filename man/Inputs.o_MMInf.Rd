% File man/Inputs.o_MMInf.Rd
\name{Inputs.o_MMInf}
\alias{Inputs.o_MMInf}
\title{Returns the input parameters of a M/M/Infinite queueing model}
\description{
  Returns the inputs parameters of a already built M/M/Infinite queueing model
}
\usage{
  \method{Inputs}{o_MMInf}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Returns the input parameters of a M/M/Infinite queueing model. The inputs parameters are created calling previously the \link{NewInput.MMInf}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MMInf}}.
}


\examples{
## create input parameters
i_mminf <- NewInput.MMInf(lambda=0.25, mu=4, n=0)

## Build the model
o_mminf <- QueueingModel(i_mminf)

## Retunns the Inputs
Inputs(o_mminf)
}

\keyword{M/M/Infinite}
