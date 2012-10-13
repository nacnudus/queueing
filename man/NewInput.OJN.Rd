% File man/NewInput.OJN.Rd
\name{NewInput.OJN}
\alias{NewInput.OJN}
\title{Define the inputs of a Open Jackson Network}
\description{
  Define the inputs of a Open Jackson Network
}

\usage{
  NewInput.OJN(prob=NULL, ...)
}
\arguments{
  \item{prob}{It is probability transition matrix or visit ratio vector. That is, the prob[i, j] is the transition probability of node i to node j, or prob[i] is the visit ratio to node i (the visit ratio values doesn't need to be probabilities, that is, a value greater than 1 can be used here. See the examples)}
  \item{\dots}{a separated by comma list of nodes of i_MM1, i_MMC or i_MMInf class}
}

\details{Define the inputs of a Open Jackson Network}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.\cr
\cr
  [Lazowska84]  Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey

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

## definition of the transition probabilities
prob <- matrix(data=c(0, 0.2, 0.56, 0.24, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE)

ojn1 <- NewInput.OJN(prob, n1, n2, n3, n4)

## Using visit ratios. Values taken from [Lazowska84], pag. 113.

## E[S] cpu = 0.005, Visit cpu = 121, D cpu = E[S] cpu * Visit cpu = 0.605
cpu <- NewInput.MM1(lambda=0.2, mu=1/0.005)

## E[S] disk1 = 0.030, Visit disk1 = 70, D disk1 = E[S] disk1 * Visit disk1 = 2.1
disk1 <- NewInput.MM1(lambda=0.2, mu=1/0.030)

## E[S] disk2 = 0.027, Visit disk2 = 50, D disk2 = E[S] disk2 * Visit disk2 = 1.35
disk2 <- NewInput.MM1(lambda=0.2, mu=1/0.027)

## In this example, to have the throughput per node, the visit ratios has to be given in this form.
## Please, don't use in the closed Jackson Network 
visit <- c(121, 70, 50)
net <- NewInput.OJN(visit, cpu, disk1, disk2)

}

\keyword{Open Jackson Network}
