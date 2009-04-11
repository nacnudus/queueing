% File man/CheckInput.i_MM1K.Rd
\name{CheckInput.i_MM1K}
\alias{CheckInput.i_MM1K}
\title{Checks the input params of a M/M/1/K queueing model}
\description{
  Checks the input params of a M/M/1/K queueing model
}
\usage{
  \method{CheckInput}{i_MM1K}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MM1K}
  \item{\dots}{aditional arguments}
}

\details{Checks the input params of a M/M/1/K queueing model. The inputs params are created calling previously the \link{NewInput.MM1K}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MM1K}}.
}

\examples{
## See example 10.7 in reference [Sixto2004] for more details.
## create input parameters
i_mm1k <- NewInput.MM1K(lambda=5, mu=5.714, k=15)

## check the parameters
CheckInput(i_mm1k)

}
\keyword{M/M/1/K}

