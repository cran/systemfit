.prepareData.systemfit <- function( data )
{
   callNoDots <- match.call( expand.dots = FALSE ) #-"- without ...-expansion

   # model frame (without formula)
   modelFrame <- callNoDots[ c( 1, match( "data", names( callNoDots ), 0 ) ) ]
   modelFrame[[1]] <- as.name( "model.frame" )

   return( modelFrame )
}




