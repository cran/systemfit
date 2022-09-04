.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0( "\nPlease cite the 'systemfit' package as:\n",
         "Arne Henningsen and Jeff D. Hamann (2007). ",
         "systemfit: A Package for Estimating Systems of Simultaneous Equations in R. ",
         "Journal of Statistical Software 23(4), 1-40. ",
         "http://www.jstatsoft.org/v23/i04/.\n\n",
         "If you have questions, suggestions, or comments ",
         "regarding the 'systemfit' package, ",
         "please use a forum or 'tracker' at systemfit's R-Forge site:\n",
         "https://r-forge.r-project.org/projects/systemfit/"),
      domain = NULL,  appendLF = TRUE )
}

# .onLoad <- function( lib, pkg ) {
#    options(Matrix.warnDeprecatedCoerce = 2) # where n >= 1
# }