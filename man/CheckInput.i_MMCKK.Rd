% File man/CheckInput.i_MMCKK.Rd
\name{CheckInput.i_MMCKK}
\alias{CheckInput.i_MMCKK}
\title{Checks the input params of a M/M/c/K/K queueing model}
\description{
  Checks the input params of a M/M/c/K/K queueing model
}
\usage{
  \method{CheckInput}{i_MMCKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMCKK}
  \item{\dots}{aditional arguments}
}

\details{Checks the input params of a M/M/c/K/K queueing model. The inputs params are created calling previously the \link{NewInput.MMCKK}}

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
i_mmckk <- NewInput.MMCKK(lambda=8, mu=2, c=5, k=12, method=0)

## check the parameters
CheckInput(i_mmckk)

}

\keyword{M/M/c/K/K}
