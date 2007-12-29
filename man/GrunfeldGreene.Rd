\name{GrunfeldGreene}
\alias{GrunfeldGreene}
\docType{data}
\title{Grunfeld Data as published by Greene (2003)}
\description{
   Panel data on 5 US firms for the years 1935-1954.
}

\usage{
   data("GrunfeldGreene")
}
\format{
A data frame containing 20 annual observations on 3 variables for 5 firms.
\describe{
   \item{invest}{gross investment.}
   \item{value}{market value of the firm (at the end of the previous year).}
   \item{capital}{capital stock of the firm (at the end of the previous year).}
   \item{firm}{name of the firm ("General Motors", "Chrysler",
      "General Electric", "Westinghouse" or "US Steel").}
   \item{year}{year.}
}
}

\details{
There exist several different versions of this data set,
and this version is considered incorrect
(see \url{http://www.stanford.edu/~clint/bench/grunfeld.htm} for details).
However, we provide this incorrect version to replicate the results
published in Theil (1971) and Greene (2003).
A correct version of this data set with 5 additional firms
is available in the \code{Ecdat} package
(data set \code{\link[Ecdat]{Grunfeld}}).
}

\source{
   Greene (2003), Appendix F, Data Sets Used in Applications, Table F13.1.
   \url{http://pages.stern.nyu.edu/~wgreene/Text/econometricanalysis.htm}
   (a subset of this data set is available in Theil (1971), p. 296).
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Prentice Hall,
Upper Saddle River (NJ).

Grunfeld, Y. (1958). \emph{The Determinants of Corporate Investment},
Unpublished Ph.D. Dissertation, University of Chicago.

Theil, Henri (1971). \emph{Principles of Econometrics},
   John Wiley & Sons, New York.
}

\examples{
## Repeating the OLS and SUR estimations in Greene (2003, pp. 351)
data( "GrunfeldGreene" )
library( plm )
GGPanel <- plm.data( GrunfeldGreene, c( "firm", "year" ) )
formulaGrunfeld <- invest ~ value + capital
# OLS
greeneOls <- systemfit( formulaGrunfeld, "OLS",
   data = GGPanel )
summary( greeneOls )
sapply( greeneOls$eq, function(x){return(summary(x)$ssr/20)} ) # sigma^2
# OLS Pooled
greeneOlsPooled <- systemfit( formulaGrunfeld, "OLS",
   data = GGPanel, pooled = TRUE )
summary( greeneOlsPooled )
sum( sapply( greeneOlsPooled$eq, function(x){return(summary(x)$ssr)}) )/97 # sigma^2
# SUR
greeneSur <- systemfit( formulaGrunfeld, "SUR",
   data = GGPanel, methodResidCov = "noDfCor" )
summary( greeneSur )
# SUR Pooled
greeneSurPooled <- systemfit( formulaGrunfeld, "SUR",
   data = GGPanel, pooled = TRUE, methodResidCov = "noDfCor",
   residCovWeighted = TRUE )
summary( greeneSurPooled )

## Repeating the OLS and SUR estimations in Theil (1971, pp. 295, 300)
GrunfeldTheil <- subset( GrunfeldGreene,
   firm \%in\% c( "General Electric", "Westinghouse" ) )
GTPanel <- plm.data( GrunfeldTheil, c( "firm", "year" ) )
formulaGrunfeld <- invest ~ value + capital
# OLS
theilOls <- systemfit( formulaGrunfeld, "OLS",
   data = GTPanel )
summary( theilOls )
# SUR
theilSur <- systemfit( formulaGrunfeld, "SUR",
   data = GTPanel, methodResidCov = "noDfCor" )
summary( theilSur )
}

\keyword{datasets}
