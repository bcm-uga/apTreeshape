\name{as.phylo.treeshape}
\alias{as.phylo.treeshape}

\title{ Conversion among tree objects }

\description{
\code{as.phylo} is a generic function - described in the APE package - which converts an object into a tree of class \code{"phylo"}. This method is an adataption of this generic method to convert objects of class \code{"treeshape"} in objects of class \code{"phylo"}.
  
  
}
\usage{
\method{as.phylo}{treeshape}(x, ...)
}

\arguments{
\item{x}{An object of class\code{"treeshape"}.}
\item{\dots}{further arguments to be passed to or from other methods.}
}

\value{
An object of class \code{"phylo"}.
}

\author{
Michael Blum <\email{michael.blum@imag.fr}>\cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}>\cr
Eric Durand <\email{eric.durand@imag.fr}>\cr
Olivier Francois <\email{olivier.francois@imag.fr}>
}

\seealso{ 
\code{\link[ape]{as.phylo}}\cr
\code{\link{as.treeshape}}\cr
}
\examples{

data(primates)
plot(primates)

library(ape)
  
primates.phylo=as.phylo(primates)
plot(primates.phylo)
}

\keyword{manip}

