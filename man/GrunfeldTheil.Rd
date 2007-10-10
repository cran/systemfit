\name{GrunfeldTheil}
\alias{GrunfeldTheil}
\docType{data}
\title{Grunfeld Data as published in Theil (1971)}
\description{
   Panel data on 2 US firms for the years 1935-1954.
}

\usage{
   data("GrunfeldTheil")
}
\format{
A data frame containing 20 annual observations on 3 variables for 2 firms.
\describe{
   \item{invest}{gross investment.}
   \item{value}{market value of the firm (at the end of the previous year).}
   \item{capital}{capital stock of the firm (at the end of the previous year).}
   \item{firm}{name of the firm ("General Electric" or "Westinghouse").}
   \item{year}{year.}
}
}

\source{
   Theil (1971), p. 296.
}

\references{
Grunfeld, Y. (1958). \emph{The Determinants of Corporate Investment},
Unpublished Ph.D. Dissertation, University of Chicago.

Theil, Henri (1971). \emph{Principles of Econometrics}, John Wiley & Sons, New York.
}

\examples{
## Repeating the OLS and SUR estimations in Theil (1971, pp. 295, 300)
data( "GrunfeldTheil" )
formulaGrunfeld <- invest ~ value + capital
# OLS
theilOls <- systemfitClassic( "OLS", formulaGrunfeld, "firm", "year",
   data = GrunfeldTheil )
summary( theilOls )
# SUR
theilSur <- systemfitClassic( "SUR", formulaGrunfeld, "firm", "year",
   data = GrunfeldTheil, rcovformula = 0 )
summary( theilSur )
}

\keyword{datasets}
