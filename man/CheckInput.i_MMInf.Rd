% File man/CheckInput.i_MMInf.Rd
\name{CheckInput.i_MMInf}
\alias{CheckInput.i_MMInf}
\title{Checks the input params of a M/M/Infinite queueing model}
\description{
  Checks the input params of a M/M/Infinite queueing model
}
\usage{
  \method{CheckInput}{i_MMInf}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Checks the input params of a M/M/Infinite queueing model. The inputs params are created calling previously the \link{NewInput.MMInf}}

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

## Check the parameters
CheckInput(i_mminf)

}
\keyword{M/M/Infinite}
