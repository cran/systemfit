
% 	$Id: ppine.Rd,v 1.2 2004/02/05 19:12:56 hamannj Exp $	


\keyword{datasets}
\name{ppine}
\alias{ppine}
\non_function{}
\title{Tree Growth Data for Ponderosa Pine}
\usage{data(ppine)}
\description{
A subset of tree growth observations from a Ponderosa pine growth database.

The \code{ppine} data frame has 166 rows and 8 columns.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{elev}{
      Altitude of the plot, in feet above mean sea level.
    }
    \item{smi}{
      Summer moisture index is the ratio of growing season heating
      degree days to growing season precipitation.
    }
    \item{dbh}{
      Diameter of the tree at breast height (4.5 feet).
    }
    \item{tht}{
      Total stem height for the tree.
    }
    \item{cr}{
      Crown ratio code. The scale is from 1 to 9 where a crown class of
      one represents a crown ratio between 0 and 15 percent. A crown
      ratio code of 2 represents a crown ratio value between 16 and
      25\%,...,8=76-85\%, 9 >=85\%.
    }
    \item{ba}{
      Plot basal area at the beginning of the growth period.
    }
    \item{dg}{
      Five-year diameter increment.
    }
    \item{hg}{
      Five-year height increment.
    }
  }
}

\details{
 The exogenous variables are \code{elev}, \code{smi}, \code{dbh},
 \code{tht}, \code{cr}, and \code{ba}; the endogenous variables
 \code{dg} and \code{hg}. There are no lagged  variables in the dataset
 and the observations are for a single remeasurement.

 The data was provided by the USDA Forest Service Intermountain
 Research Station.

 
%  The data was provided by the USDA Forest Service Intermountain
%  Research Station from installations on USDA Forest Service, Spokane Indian
%  Reservation and the Nez Perce Indian Reservation lands.
}

\source{
  William R. Wykoff \email{wwykoff@fs.fed.us} \emph{Rocky Mountain Research Station, 1221 South Main Street, Moscow, ID 83843}
}


\examples{
data(ppine)
}

\keyword{datasets}


