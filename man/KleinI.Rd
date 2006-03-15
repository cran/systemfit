\name{KleinI}
\alias{KleinI}
\title{Klein Model I}
\description{
   Data for Klein's (1950) Model I of the US economy.
}

\usage{data("KleinI")}

\format{
A data frame containing annual observations from 1920 to 1941
   \describe{
      \item{year}{Year.}
      \item{consump}{Consumption.}
      \item{corpProf}{Corporate profits.}
      \item{corpProfLag}{Corporate profits of the previous year.}
      \item{privWage}{Private wage bill.}
      \item{invest}{Investment.}
      \item{capitalLag}{Capital stock of the previous year.}
      \item{gnp}{Gross national product.}
      \item{gnpLag}{Gross national product of the previous year.}
      \item{govWage}{Government wage bill.}
      \item{govExp}{Government spending.}
      \item{taxes}{Taxes.}
      \item{wages}{Sum of private and government wage bill.}
      \item{trend}{time trend measured as years from 1931.}
   }
}

\source{
   Greene (2003), Appendix F, Data Sets Used in Applications, Table F15.1.

   \url{http://pages.stern.nyu.edu/~wgreene/Text/econometricanalysis.htm}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Prentice Hall,
Upper Saddle River (NJ).

Klein, L. (1950). \emph{Economic Fluctuations in the United States, 1921--1941}.
John Wiley, New York.
}

\examples{
## Repeating the estimations of Klein's (1950) Model I
## in Greene (2003, pp. 381 and 412)
data( "KleinI" )
eqConsump  <- consump ~ corpProf + corpProfLag + wages
eqInvest   <- invest ~ corpProf + corpProfLag + capitalLag
eqPrivWage <- privWage ~ gnp + gnpLag + trend
inst <- ~ govExp + taxes + govWage + trend + capitalLag + corpProfLag + gnpLag
system <- list( Consumption = eqConsump, Investment = eqInvest,
   "Private Wages" = eqPrivWage )
# OLS
kleinOls <- systemfit( "OLS", system, data = KleinI )
summary( kleinOls )
# 2SLS
klein2sls <- systemfit( "2SLS", system, inst = inst, data = KleinI,
   rcovformula = 0 )
summary( klein2sls )
# 3SLS
klein3sls <- systemfit( "3SLS", system, inst = inst, data = KleinI,
   rcovformula = 0 )
summary( klein3sls )
# I3SLS
kleinI3sls <- systemfit( "3SLS", system, inst = inst, data = KleinI,
   rcovformula = 0, maxit = 500 )
summary( kleinI3sls )
}
\keyword{datasets}
