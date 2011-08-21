% File man/summary.o_MCCN.Rd
\name{summary.o_MCCN}
\alias{summary.o_MCCN}
\title{Reports the results of a MultiClass Closed Network}
\description{
  Reports the results of a MultiClass Closed Network
}
\usage{
  \method{summary}{o_MCCN}(object, \dots)
}
\arguments{
  \item{object}{a object of class o_MCCN}
  \item{\dots}{aditional arguments}
}

\details{generates a report of the queueing network received as parameter}

\references{
  [Lazowska84] Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey
}

\seealso{
  \code{\link{QueueingModel.i_MCCN}}.
}


\examples{
## See example in pag 142 in reference [Lazowska84] for more details.

classes <- 2
vNumber <- c(1, 1)
vThink <- c(0, 0)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=c(10, 9, 5, 4), nrow=2, ncol=2, byrow=TRUE)
vMu <- matrix(data=c(10, 3, 5/2, 1), nrow=2, ncol=2, byrow=TRUE)

i_MCCN1 <- NewInput.MCCN(classes, vNumber, vThink, nodes, vType, vVisit, vMu)

# Build the model
o_MCCN1 <- QueueingModel(i_MCCN1)

summary(o_MCCN1)

}
\keyword{MultiClass Closed Network}
