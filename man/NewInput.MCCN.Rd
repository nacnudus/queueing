% File man/NewInput.MCCN.Rd
\name{NewInput.MCCN}
\alias{NewInput.MCCN}
\title{Define the inputs of a MultiClass Closed Network}
\description{
  Define the inputs of a MultiClass Closed Network
}

\usage{
NewInput.MCCN(
  classes, vNumber, vThink, nodes, vType, vVisit, vService, method=1, tol=0.01
)
}

\arguments{
  \item{classes}{The number of classes}
  \item{vNumber}{A vector with the number of customers of each class}
  \item{vThink}{A vector with the think time of each class}
  \item{nodes}{The number of nodes in the network}
  \item{vType}{A vector with the type of node: "Q" for queueing nodes or "D" for delay nodes}
  \item{vVisit}{A matrix[i, j]. The rows represents the different visit count for each class i to each node j}
  \item{vService}{A matrix[i, j]. The rows represents the different service time for each class i in each node j}
  \item{method}{If method is 0, the exact MVA algorith is used. If method is 1, the Bard-Schweitzer approximation algorithm is used}
  \item{tol}{If the parameter method is 1, this is the tolerance parameter of the algorithm}
}

\details{Define the inputs of a MultiClass Closed Network}

\references{
  [Lazowska84]  Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey

}

\seealso{
  \code{\link{QueueingModel.i_MCCN}}
}


\examples{
## See example in pag 142 in reference [Lazowska84] for more details.

classes <- 2
vNumber <- c(1, 1)
vThink <- c(0, 0)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=c(10, 9, 5, 4), nrow=2, ncol=2, byrow=TRUE)
vService <- matrix(data=c(1/10, 1/3, 2/5, 1), nrow=2, ncol=2, byrow=TRUE)

i_MCCN1 <- NewInput.MCCN(classes, vNumber, vThink, nodes, vType, vVisit, vService)

}

\keyword{MultiClass Closed Network}
