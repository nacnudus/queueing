\name{queueing-package}
\alias{queueing-package}
\alias{queueing}
\docType{package}

\title{
Analysis of Queueing Networks and Models.
}

\description{
It provides a versatile tool for analysis of birth and death based Markovian Queueing Models and Single and Multiclass Product-Form Queueing Networks.\cr\cr
It implements the following basic markovian models:\cr\cr

\tabular{llll}{
M/M/1, \tab M/M/c, \tab \tab M/M/Infinite,\cr
M/M/1/K, \tab M/M/c/K, \tab M/M/c/c, \tab\cr
M/M/1/K/K, \tab M/M/c/K/K, \tab M/M/c/K/m, \tab M/M/Infinite/K/K\cr\cr
}


It also solves the following types of networks:\cr

\itemize{
\item Multiple Channel Open Jackson Networks.
\item Multiple Channel Closed Jackson Networks.
\item Single Channel Multiple Class Open Networks.
\item Single Channel Multiple Class Closed Networks
\item Single Channel Multiple Class Mixed Networks 
}

Also it provides B-Erlang, C-Erlang and Engset calculators.\cr

This work is dedicated to the memory of D. Sixto Rios Insua.
}

\details{

All models are used in the same way:\cr

\enumerate{
\item Create inputs calling the appropiate \emph{NewInput.model}. For example, \code{x <- NewInput.MM1(lambda=0.25, mu=1, n=10)} for a M/M/1 model. To know the exact acronymn model to use for \emph{NewInput} function, you can search the html help or write \code{help.search("NewInput")} at the command line.
\item Optionally, as a help for creating the inputs, the \code{CheckInput(x)} function can be called
\item Solve the model calling \code{y <- QueueingModel(x)}. In this step, the \code{CheckInput(x)} will be called. That is the reason that the previous step is optional
\item Finally, you can get a performance value as \code{W(y)}, \code{Wq(y)} or a report of the principals performace values calling \code{summary(y)}
}


See the examples for more detailed information of the use.

}

\author{
Author, Maintainer and Copyright: Pedro Canadilla \email{pedro.canadilla@gmail.com}
}

\references{
[Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\keyword{ queueing }

\examples{

## M/M/1 model
summary(QueueingModel(NewInput.MM1(lambda=1/4, mu=1/3, n=0)))

## M/M/1/K model
summary(QueueingModel(NewInput.MM1K(lambda=1/4, mu=1/3, k=3)))

}
