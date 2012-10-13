% File man/Throughput.o_MCON.Rd
\name{Throughput.o_MCON}
\alias{Throughput.o_MCON}
\title{Reports the throughput of a MultiClass Open Network}
\description{
  Reports the throughput of a MultiClass Open Network
}
\usage{
  \method{Throughput}{o_MCON}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MCON}
  \item{\dots}{aditional arguments}
}

\details{
Reports the throughput of a MultiClass Open Network
}

\references{
  [Lazowska84] Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey
}

\seealso{
  \code{\link{NewInput.MCON}}, \code{\link{CheckInput.i_MCON}}, \code{\link{QueueingModel.i_MCON}}
}


\examples{
## See example in pag 138 in reference [Lazowska84] for more details.

classes <- 2
vLambda <- c(3/19, 2/19)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=c(10, 9, 5, 4), nrow=2, ncol=2, byrow=TRUE)
vService <- matrix(data=c(1/10, 1/3, 2/5, 1), nrow=2, ncol=2, byrow=TRUE)

i_mcon1 <- NewInput.MCON(classes, vLambda, nodes, vType, vVisit, vService)

# Build the model
o_mcon1 <- QueueingModel(i_mcon1)

Throughput(o_mcon1)

}
\keyword{MultiClass Open Network}

