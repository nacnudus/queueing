% File man/NewInput.MCON.Rd
\name{NewInput.MCON}
\alias{NewInput.MCON}
\title{Define the inputs of a new one class of MultiClass Open Network}
\description{
  Define the inputs of a new one class of MultiClass Open Network
}

\usage{
  \method{NewInput}{MCON}(vLambda, nodes, vType, vVisit, vMu)
}
\arguments{
  \item{vLambda}{It is a vector with the rate of arrivals of each class}
  \item{nodes}{the number of nodes in the network}
  \item{vType}{A vector with the type of node: "Q" for queueing nodes or "D" for delay nodes}
  \item{vVisit}{A matrix[i, j]. The rows represents the different visit count for each class i to each node j}
  \item{vMu}{A matrix[i, j]. The rows represents the different service rate for each class i in each node j}  
}

\details{Define the inputs of a new one class MultiClass Open Network}

\references{
  [Lazowska84]  Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey

}

\seealso{
  \code{\link{QueueingModel.i_MCON}}
}


\examples{
## See example in pag 138 in reference [Lazowska84] for more details.

vLambda <- c(3/19, 2/19)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=c(10, 9, 5, 4), nrow=2, ncol=2, byrow=TRUE)
vMu <- matrix(data=c(10, 3, 5/2, 1), nrow=2, ncol=2, byrow=TRUE)

i_mcon1 <- NewInput.MCON(vLambda, nodes, vType, vVisit, vMu)

}

\keyword{MultiClass Open Network}
