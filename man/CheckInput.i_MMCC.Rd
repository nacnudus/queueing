% File man/CheckInput.i_MMCC.Rd
\name{CheckInput.i_MMCC}
\alias{CheckInput.i_MMCC}
\title{Checks the input params of a M/M/c/c queueing model}
\description{
  Checks the input params of a M/M/c/c queueing model
}
\usage{
  \method{CheckInput}{i_MMCC}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMCC}
  \item{\dots}{aditional arguments}
}

\details{Checks the input params of a M/M/c/c queueing model. The inputs params are created calling previously the \link{NewInput.MMCC}}

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

## check the parameters
CheckInput(i_mmcc)

}
\keyword{M/M/c/c}

