% File man/CheckInput.i_MMCKM.Rd
\name{CheckInput.i_MMCKM}
\alias{CheckInput.i_MMCKM}
\title{Checks the input params of a M/M/c/K/m queueing model}
\description{
  Checks the input params of a M/M/c/K/m queueing model
}
\usage{
  \method{CheckInput}{i_MMCKM}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMCKM}
  \item{\dots}{aditional arguments}
}

\details{Checks the input params of a M/M/c/K/m queueing model. The inputs params are created calling previously the \link{NewInput.MMCKM}}

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
i_mmckm <- NewInput.MMCKM(lambda=0.25, mu=4, c=2, k=4, m=8, method=0)

## check the parameters
CheckInput(i_mmckm)

}
\keyword{M/M/c/K/m}
