% File man/summary.o_MCON.Rd
\name{summary.o_MCON}
\alias{summary.o_MCON}
\title{Reports the results of a MultiClass Open Network}
\description{
  Reports the results of a MultiClass Open Network
}
\usage{
  \method{summary}{o_MCON}(object, \dots)
}
\arguments{
  \item{object}{a object of class o_MCON}
  \item{\dots}{aditional arguments}
}

\details{generates a report of the queueing network received as parameter}

\references{
  [Lazowska84] Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey
}

\seealso{
  \code{\link{QueueingModel.i_MCON}}.
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

summary(o_mcon1)

}
\keyword{MultiClass Open Network}
