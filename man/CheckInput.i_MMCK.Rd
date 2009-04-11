% File man/CheckInput.i_MMCK.Rd
\name{CheckInput.i_MMCK}
\alias{CheckInput.i_MMCK}
\title{Checks the input params of a M/M/c/K queueing model}
\description{
  Checks the input params of a M/M/c/K queueing model
}
\usage{
  \method{CheckInput}{i_MMCK}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMCK}
  \item{\dots}{aditional arguments}
}

\details{Checks the input params of a M/M/c/K queueing model. The inputs params are created calling previously the \link{NewInput.MMCK}}

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

## Check the inputs
CheckInput(i_mmck)

}
\keyword{M/M/c/K}

