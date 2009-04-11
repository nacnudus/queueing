% File man/CheckInput.Rd
\name{CheckInput}
\alias{CheckInput}
\title{Generic S3 method to check the params of a queueing model (or network)}
\description{
  Generic S3 method to check the params of a queueing model (or network)
}
\usage{
  CheckInput(x, \dots)
}
\arguments{
  \item{x}{a object of class i_MM1, i_MMC, i_MM1K, i_MMCK, i_MM1KK, i_MMCKK, i_MMCC, i_MMCKM, i_MMInf, i_OJN}
  \item{\dots}{aditional arguments}
}

\details{Generic S3 method to check the params of a queueing model (or network)}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MM1}}\cr
  \code{\link{CheckInput.i_MMC}}\cr
  \code{\link{CheckInput.i_MM1K}}\cr
  \code{\link{CheckInput.i_MMCK}}\cr
  \code{\link{CheckInput.i_MM1KK}}\cr
  \code{\link{CheckInput.i_MMCKK}}\cr
  \code{\link{CheckInput.i_MMCC}}\cr
  \code{\link{CheckInput.i_MMCKM}}\cr
  \code{\link{CheckInput.i_MMInf}}\cr
  \code{\link{CheckInput.i_OJN}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Check the inputs
CheckInput(i_mm1)

}
\keyword{Queueing Models}
