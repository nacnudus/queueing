% File man/NewInput.MMCC.Rd
\name{NewInput.MMCC}
\alias{NewInput.MMCC}
\title{Define the inputs of a new M/M/c/c queueing model}
\description{
  Define the inputs of a new M/M/c/c queueing model
}

\usage{
NewInput.MMCC(lambda=0, mu=0, c=1, method=1)
}

\arguments{
  \item{lambda}{arrival rate}
  \item{mu}{server service rate}
  \item{c}{number of servers}
  \item{method}{with method = 0, the state probabilities are calculated using the formal definition (with overflow problems with factorials; with method = 1 (default), the truncated poisson distribution is used (recomended for professional use)}
}

\details{Define the inputs of a new M/M/c/c queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.\cr
\cr
  [Kobayashi2012] Hisashi Kobayashi, Brian L. Mark, William Turin (2012).\cr
\emph{Probability, Random Processes, and Statistical Analysis: Applications to Communications, Signal Processing, Queueing Theory and Mathematical Finance}.\cr
Cambridge University Press.
}

\seealso{
  \code{\link{CheckInput.i_MMCC}}
}



\examples{
## See example 10.12 in reference [Sixto2004] for more details.
## create input parameters
i_mmcc <- NewInput.MMCC(lambda=3, mu=0.25, c=15)

}
\keyword{M/M/c/c}

