% File man/CheckInput.i_MCON.Rd
\name{CheckInput.i_MCON}
\alias{CheckInput.i_MCON}
\title{Check the input params of a MultiClass Open Network}
\description{
  Check the input params of a MultiClass Open Network
}
\usage{
  \method{CheckInput}{i_MCON}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MCON}
  \item{\dots}{aditional arguments}
}

\details{Check the input params of a MultiClass Open Network. The inputs params are created calling previously the \link{NewInput.MCON}}

\references{
  [Lazowska84] Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey
}

\seealso{
  \code{\link{NewInput.MCON}}
}


\examples{
## See example in pag 138 in reference [Lazowska84] for more details.

vLambda <- c(3/19, 2/19)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=c(10, 9, 5, 4), nrow=2, ncol=2, byrow=TRUE)
vMu <- matrix(data=c(10, 3, 5/2, 1), nrow=2, ncol=2, byrow=TRUE)

o_mcon1 <- NewInput.MCON(vLambda, nodes, vType, vVisit, vMu)

CheckInput(o_mcon1)

}
\keyword{MultiClass Open Network}
