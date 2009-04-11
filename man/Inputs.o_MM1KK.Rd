% File man/Inputs.o_MM1KK.Rd
\name{Inputs.o_MM1KK}
\alias{Inputs.o_MM1KK}
\title{Returns the input parameters of a M/M/1/K/K queueing model}
\description{
  Returns the inputs parameters of a already built M/M/1/K/K queueing model
}
\usage{
  \method{Inputs}{o_MM1KK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{Returns the input parameters of a M/M/1/K/K queueing model. The inputs parameters are created calling previously the \link{NewInput.MM1KK}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MM1KK}}.
}


\examples{
## See example 10.13 in reference [Sixto2004] for more details.
## create input parameters
i_mm1kk <- NewInput.MM1KK(lambda=0.25, mu=4, k=2)

## Build the model
o_mm1kk <- QueueingModel(i_mm1kk)

## Retunns the Inputs
Inputs(o_mm1kk)
}

\keyword{M/M/1/K/K}
