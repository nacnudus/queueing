% File man/Wc.o_MCCN.Rd
\name{Wc.o_MCCN}
\alias{Wc.o_MCCN}
\title{Returns the vector with each class mean time spend on a MultiClass Closed Network}
\description{
  Returns the vector with each class mean time spend on a MultiClass Closed Network
}

\usage{
  \method{Wc}{o_MCCN}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MCCN}
  \item{\dots}{aditional arguments}
}

\details{Returns the vector with each class mean time spend on a MultiClass Closed Network}

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

Wc(o_MCCN1)

}

\keyword{MultiClass Closed Network}
