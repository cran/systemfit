\name{GrunfeldGreene}
\alias{GrunfeldGreene}
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
However, we provide this incorrect version to replicated the results
published in Greene (2003).
}

\source{
   Greene (2003), Appendix F, Data Sets Used in Applications, Table F13.1.

   \url{http://pages.stern.nyu.edu/~wgreene/Text/econometricanalysis.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Prentice Hall,
Upper Saddle River (NJ).

Grunfeld, Y. (1958). \emph{The Determinants of Corporate Investment},
Unpublished Ph.D. Dissertation, University of Chicago.
}

\examples{
## Repeating the OLS and SUR estimations in Greene (2003, pp. 351)
data( "GrunfeldGreene" )
formulaGrunfeld <- invest ~ value + capital
# OLS
greeneOls <- systemfitClassic( "OLS", formulaGrunfeld, "firm", "year",
   data = GrunfeldGreene )
summary( greeneOls )
sapply( greeneOls$eq, function(x){return(x$ssr/20)} ) # sigma^2
# OLS Pooled
greeneOlsPooled <- systemfitClassic( "OLS", formulaGrunfeld, "firm", "year",
   data = GrunfeldGreene, pooled = TRUE )
summary( greeneOlsPooled )
sum( sapply( greeneOlsPooled$eq, function(x){return(x$ssr)}) )/97 # sigma^2
# SUR
greeneSur <- systemfitClassic( "SUR", formulaGrunfeld, "firm", "year",
   data = GrunfeldGreene, rcovformula = 0 )
summary( greeneSur )
# SUR Pooled
greeneSurPooled <- systemfitClassic( "WSUR", formulaGrunfeld, "firm", "year",
   data = GrunfeldGreene, pooled = TRUE, rcovformula = 0 )
summary( greeneSurPooled )
}

\keyword{datasets}
