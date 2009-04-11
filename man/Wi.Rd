% File man/Wi.Rd
\name{Wi}
\alias{Wi}
\title{Generic S3 method to return the mean time spend in each node (or server) of a network}
\description{
  Generic S3 method to return the mean time spend in each node (or server) of a network
}
\usage{
  Wi(x, \dots)
}
\arguments{
  \item{x}{a object of class o_OJN, o_CJN}
  \item{\dots}{aditional arguments}
}

\details{Generic S3 method to return the mean time spend in each node (or server) of a network}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{Wi.o_OJN}}\cr
  \code{\link{Wi.o_CJN}}
}


\examples{
## See example 11.11 in reference [Sixto2004] for more details.
## create the nodes
n1 <- NewInput.MM1(lambda=8, mu=14, n=0)
n2 <- NewInput.MM1(lambda=0, mu=8, n=0)
n3 <- NewInput.MM1(lambda=6, mu=17, n=0)
n4 <- NewInput.MM1(lambda=0, mu=7, n=0)

# definition of the transition probabilities
prob <- matrix(data=c(0, 0.2, 0.56, 0.24, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)

ojn1 <- NewInput.OJN(prob, n1, n2, n3, n4)

m_ojn1 <- QueueingModel(ojn1)

Wi(m_ojn1)


}
\keyword{Queueing Models}
