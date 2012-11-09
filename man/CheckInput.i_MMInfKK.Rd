% File man/CheckInput.i_MMInfKK.Rd
\name{CheckInput.i_MMInfKK}
\alias{CheckInput.i_MMInfKK}
\title{Checks the input params of a M/M/Infinite/K/K queueing model}
\description{
  Checks the input params of a M/M/Infinite/K/K queueing model
}
\usage{
  \method{CheckInput}{i_MMInfKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MMInfKK}
  \item{\dots}{aditional arguments}
}

\details{Checks the input params of a M/M/Infinite/K/K queueing model. The inputs params are created calling previously the \link{NewInput.MMInfKK}}

\references{
  [Kleinrock1975] Leonard Kleinrock (1975).\cr
\emph{Queueing Systems Vol 1: Theory}.\cr
John Wiley & Sons.
}

\seealso{
  \code{\link{NewInput.MMInfKK}}.
}


\examples{
## create input parameters
i_MMInfKK <- NewInput.MMInfKK(lambda=0.25, mu=4, k=4)

## check the parameters
CheckInput(i_MMInfKK)

}
\keyword{M/M/Infinite/K/K}
