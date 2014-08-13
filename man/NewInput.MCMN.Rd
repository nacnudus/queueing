% File man/NewInput.MCMN.Rd
\name{NewInput.MCMN}
\alias{NewInput.MCMN}
\title{Define the inputs of a MultiClass Mixed Network}
\description{
  Define the inputs of a MultiClass Mixed Network
}

\usage{
NewInput.MCMN(
  classes, vLambda, vNumber, vThink, nodes, vType, vVisit, vService, method=0, tol=0.01
)
}

\arguments{
  \item{classes}{The number of classes}
  \item{vLambda}{It is a vector with the rate of arrivals of each class}
  \item{vNumber}{A vector with the number of customers of each class}
  \item{vThink}{A vector with the think time of each class}
  \item{nodes}{The number of nodes in the network}
  \item{vType}{A vector with the type of node: "Q" for queueing nodes or "D" for delay nodes}
  \item{vVisit}{A matrix[i, j]. The rows represents the different visit count for each class i to each node j. Take caution about the orden: open classes are defined first and closed clasess are defined second}
  \item{vService}{A matrix[i, j]. The rows represents the different service times for each class i in each node j.  Take caution about the orden: open classes are defined first and closed clasess are defined second.}
  \item{method}{If method is 0, the exact MVA algorith is used. If method is 1, the Bard-Schweitzer approximation algorithm is used}
  \item{tol}{If the parameter method is 1, this is the tolerance parameter of the algorithm}
}

\details{Define the inputs of a MultiClass Mixed Network}

\references{
  [Lazowska84]  Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey

}

\seealso{
  \code{\link{QueueingModel.i_MCMN}}
}


\examples{
## See example in pag 147 in reference [Lazowska84] for more details.

classes <- 4 # A and B are open classes and C and D are closed classes.
vLambda <- c(1, 1/2)
vNumber <- c(1, 1)
vThink <- c(0, 0)
nodes <- 2
vType <- c("Q", "Q")

# When the visit ratios and vService are set,
#   be sure that the open classes are in the first positions
#   and the closed classes after the open classes.
vVisit <- matrix(data=1, nrow=4, ncol=2)

# A and B are open clasess:
#   with demand service of 1/4 and 1/2 at the node 1 and 1/2 and 1 at the node 2
# C and D are open clasess:
#   with demand service of 1/4 and 1/2 at the node 1 and 1/2 and 1 at the node 2 
vService <- matrix(data=c(1/4, 1/2, 1/2, 1, 1/6, 1, 1, 4/3), nrow=4, ncol=2)

i_mcmn1 <- NewInput.MCMN(classes, vLambda, vNumber, vThink, nodes, vType, vVisit, vService)

}

\keyword{MultiClass Mixed Network}
