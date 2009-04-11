% File man/QueueingModel.i_MCON.Rd
\name{QueueingModel.i_MCON}
\alias{QueueingModel.i_MCON}
\title{Builds one class MultiClass Open Network}
\description{
  Builds one class MultiClass Open Network
}
\usage{
  \method{QueueingModel}{i_MCON}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MCON}
  \item{\dots}{aditional arguments}
}

\details{Build one class MultiClass Open Network. It also checks the input params calling the \link{CheckInput.i_MCON}}

\references{
  [Lazowska84] Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey
}

\seealso{
  \code{\link{CheckInput.i_MCON}}
}


\examples{
## See example in pag 138 in reference [Lazowska84] for more details.

vLambda <- c(3/19, 2/19)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=c(10, 9, 5, 4), nrow=2, ncol=2, byrow=TRUE)
vMu <- matrix(data=c(10, 3, 5/2, 1), nrow=2, ncol=2, byrow=TRUE)

i_mcon1 <- NewInput.MCON(vLambda, nodes, vType, vVisit, vMu)

# Build the model
o_mcon1 <- QueueingModel(i_mcon1)

o_mcon1

}
\keyword{MultiClass Open Network}
