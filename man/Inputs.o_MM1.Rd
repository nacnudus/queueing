% File man/Inputs.o_MM1.Rd
\name{Inputs.o_MM1}
\alias{Inputs.o_MM1}
\title{Returns the input parameters of a M/M/1 queueing model}
\description{
  Returns the inputs parameters of a already built M/M/1 queueing model
}
\usage{
  \method{Inputs}{o_MM1}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1}
  \item{\dots}{aditional arguments}
}

\details{Returns the input parameters of a M/M/1 queueing model. The inputs parameters are created calling previously the \link{NewInput.MM1}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MM1}}.
}


\examples{
## See example 10.3 in reference [Sixto2004] for more details.
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Retunns the Inputs
Inputs(o_mm1)
}

\keyword{M/M/1}
