
% 	$Id: systemfit-internal.Rd,v 1.2 2004/02/04 02:34:28 hamannj Exp $	

\name{systemfit-internal}
\alias{knls}

\title{Internal systemfit functions}
\description{
  Internal systemfit functions
}
\usage{
knls( theta, eqns, data, fitmethod="OLS", parmnames, instr=NULL, S=NULL )
}

\details{
  These functions are no the called by the user.
  The knls is the driver function that is passed to the nlm call in
  \code{\link{nlsystemfit}}. 
}

%\seealso{\code{\link{nlm}},\code{\link{qr}}, and \code{\link{systemfit}} }


\keyword{internal}
