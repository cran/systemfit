\keyword{datasets}
\name{kmenta}
\alias{kmenta}
\non_function{}
\title{Partly Artificial Data on the U. S. Economy}
\usage{data(kmenta)}
\description{
These are partly contrived data from Kmenta (1986), constructed
to illustrate estimation of a simultaneous-equation model.

The \code{kmenta} data frame has 20 rows and 5 columns.
}
\format{
This data frame contains the following columns:
\describe{
\item{q}{
food consumption per capita.
}
\item{p}{
ratio of food prices to general consumer prices.
}
\item{d}{
disposable income in constant dollars.
}
\item{f}{
ratio of preceding year's prices received by farmers
to general consumer prices.
}
\item{a}{
time in years.
}
}
}

\details{
 The exogenous variables \code{d}, \code{f}, and \code{a} are based on
 real data; the endogenous variables \code{p} and \code{q} were generated
 by simulation.
}
\source{
 Kmenta, J. (1986)
 \emph{Elements of Econometrics, Second Edition}, Macmillan.
}
\examples{
data(kmenta)
}
\keyword{datasets}


