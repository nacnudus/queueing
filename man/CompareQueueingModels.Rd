% File man/CompareQueueingModels.Rd
\name{CompareQueueingModels}
\alias{CompareQueueingModels}
\alias{CompareQueueingModels2}
\title{Compare several queueing models in a tabulated format}
\description{
  Compare several queueing models in a tabulated format
}
\usage{
  CompareQueueingModels(model, \dots)
  CompareQueueingModels2(models)
}
\arguments{
  \item{model}{A Queueing Model obtained calling QueueingModel from classes described in the details section}
  \item{\dots}{a separated by comma list of queueing models obtained calling QueueingModel from classes described in the details section}
  \item{models}{A list of queueing models obtained calling QueueingModel from classes described in the details section}
}

\details{Compare several queueing models in a tabulated format. By now, only o_MM1, o_MMC, o_MMInf, o_MM1K, o_MMCK, o_MMCC, o_MM1KK, o_MMCKK, o_MMCKM, o_MMInfKK classes can be compared}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel}}
}


\examples{

q1 <- QueueingModel(NewInput.MM1(lambda=5, mu=7))
q2 <- QueueingModel(NewInput.MMC(lambda=5, mu=3, c=4))
q3 <- QueueingModel(NewInput.MMInf(lambda=3, mu=4))
q4 <- QueueingModel(NewInput.MMCC(lambda=5, mu=3, c=4))

CompareQueueingModels(q1, q2, q3)
CompareQueueingModels2(list(q1, q2, q3, q4))

}
\keyword{CompareQueueingModels}
