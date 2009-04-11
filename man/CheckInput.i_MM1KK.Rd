% File man/CheckInput.i_MM1KK.Rd
\name{CheckInput.i_MM1KK}
\alias{CheckInput.i_MM1KK}
\title{Checks the input params of a M/M/1/K/K queueing model}
\description{
  Checks the input params of a M/M/1/K/K queueing model
}
\usage{
  \method{CheckInput}{i_MM1KK}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MM1KK}
  \item{\dots}{aditional arguments}
}

\details{Checks the input params of a M/M/1/K/K queueing model. The inputs params are created calling previously the \link{NewInput.MM1KK}}

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

## check the parameters
CheckInput(i_mm1kk)

}
\keyword{M/M/1/K/K}
