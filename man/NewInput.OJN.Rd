% File man/NewInput.OJN.Rd
\name{NewInput.OJN}
\alias{NewInput.OJN}
\title{Define the inputs of a new one class Open Jackson Network}
\description{
  Define the inputs of a new one class Open Jackson Network
}

\usage{
  \method{NewInput}{OJN}(prob=NULL, ...)
}
\arguments{
  \item{prob}{It is probability transition matrix or visit ratio vector. That is, the prob[i, j] is the transition probability of node i to node j, or prob[i] is the visit ratio to node i}
  \item{\dots}{a separated by comma list of nodes of i_MM1 or i_MMC class}
}

\details{Define the inputs of a new one class Open Jackson Network}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_OJN}}
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

ojn1 <- NewInput.OJN(prob, n1, n2, n3, n4)

}

\keyword{Open Jackson Network}
