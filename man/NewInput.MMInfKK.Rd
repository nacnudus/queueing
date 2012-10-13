% File man/NewInput.MMInfKK.Rd
\name{NewInput.MMInfKK}
\alias{NewInput.MMInfKK}
\title{Define the inputs of a new M/M/Inf/K/K queueing model}
\description{
  Define the inputs of a new M/M/Inf/K/K queueing model
}

\usage{
  NewInput.MMInfKK(lambda=0, mu=0, k=1)
}
\arguments{
  \item{lambda}{arrival rate}
  \item{mu}{server service rate}
  \item{k}{system capacity}
}

\details{Define the inputs of a new M/M/Inf/K/K queueing model}

\references{
  [Kleinrock1975] Leonard Kleinrock (1975).\cr
\emph{Queueing Systems Vol 1: Theory}.\cr
John Wiley & Sons.
}

\seealso{
  \code{\link{CheckInput.i_MMInfKK}}
}


\examples{
## create input parameters
i_MMInfKK <- NewInput.MMInfKK(lambda=0.25, mu=4, k=4)

}
\keyword{M/M/Inf/K/K}
