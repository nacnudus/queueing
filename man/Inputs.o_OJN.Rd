% File man/Inputs.o_OJN.Rd
\name{Inputs.o_OJN}
\alias{Inputs.o_OJN}
\title{Returns the input params of an Open Jackson Network}
\description{
  Returns the input params of an Open Jackson Network
}
\usage{
  \method{Inputs}{o_OJN}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_OJN}
  \item{\dots}{aditional arguments}
}

\details{Returns the input params of an Open Jackson Network. The inputs parameters are created calling previously the \link{NewInput.OJN}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.OJN}}.
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

i_ojn1 <- NewInput.OJN(prob, n1, n2, n3, n4)

# Build the model
o_ojn1 <- QueueingModel(i_ojn1)

Inputs(o_ojn1)

}

\keyword{Open Jackson Network}
