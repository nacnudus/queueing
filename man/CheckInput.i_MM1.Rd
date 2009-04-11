% File man/CheckInput.i_MM1.Rd
\name{CheckInput.i_MM1}
\alias{CheckInput.i_MM1}
\title{Checks the input params of a M/M/1 queueing model}
\description{
  Checks the input params of a M/M/1 queueing model
}
\usage{
  \method{CheckInput}{i_MM1}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MM1}
  \item{\dots}{aditional arguments}
}

\details{Checks the input params of a M/M/1 queueing model. The inputs params are created calling previously the \link{NewInput.MM1}}

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

## Check the inputs
CheckInput(i_mm1)

}
\keyword{M/M/1}
