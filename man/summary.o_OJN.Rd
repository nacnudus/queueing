% File man/summary.o_OJN.Rd
\name{summary.o_OJN}
\alias{summary.o_OJN}
\title{Reports the results of an Open Jackson Network}
\description{
  Reports the results of an Open Jackson Network
}
\usage{
  \method{summary}{o_OJN}(object, \dots)
}
\arguments{
  \item{object}{a object of class o_OJN}
  \item{\dots}{aditional arguments}
}

\details{generates a report of the queueing network received as parameter}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_OJN}}.
}


\examples{
## See example 11.11 in reference [Sixto2004] for more details.
## create the nodes
n1 <- NewInput.MM1(lambda=8, mu=14, n=0)
n2 <- NewInput.MM1(lambda=0, mu=9, n=0)
n3 <- NewInput.MM1(lambda=6, mu=17, n=0)
n4 <- NewInput.MM1(lambda=0, mu=7, n=0)

# definition of the transition probabilities
prob <- matrix(data=c(0, 0.2, 0.56, 0.24, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)

i_ojn <- NewInput.OJN(prob, n1, n2, n3, n4)

o_ojn <- QueueingModel(i_ojn)

summary(o_ojn)

}
\keyword{Open Jackson Network}
