% File man/CheckInput.i_MMC.Rd
\name{CheckInput.i_MMC}
\alias{CheckInput.i_MMC}
\title{Checks the input params of a M/M/c queueing model}
\description{
  Checks the input params of a M/M/c queueing model
}
\usage{
  \method{CheckInput}{i_MMC}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMC}
  \item{\dots}{aditional arguments}
}

\details{Checks the input params of a M/M/c queueing model. The inputs params are created calling previously the \link{NewInput.MMC}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MMC}}.
}


\examples{
## See example 10.9 in reference [Sixto2004] for more details.
## create input parameters
i_mmc <- NewInput.MMC(lambda=5, mu=10, c=2, n=0, method=0)

## check the parameters
CheckInput(i_mmc)

}
\keyword{M/M/c}

