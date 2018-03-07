\name{plot.treeshape}
\alias{plot.treeshape}

\title{ Plot phylogenetic treeshapes.  }

\description{
Plot method for objects of class \code{"treeshape"}.
}

\usage{
\method{plot}{treeshape}(x, y, \dots)
}

\arguments{
\item{x}{An object of class \code{"treeshape"}.}
\item{y}{An object of class \code{"treeshape"}.}
\item{...}{Further arguments to be passed to \code{plot()}.}
}

\value{
A null value is returned. Results are displayed on graphical window.
}

\details{ 
If two trees are specified, they are plotted on the same window. This option is provided is order to facilitate the comparison between two trees.
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}>\cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}>\cr
Eric Durand <\email{eric.durand@imag.fr}>\cr
Olivier Francois <\email{olivier.francois@imag.fr}>
}
 
\seealso{ 
\code{\link[graphics]{plot}} for the basic plotting function in R 
}

\examples{

## Visual representation of the universal tree of life provided in data
data(universal.treeshape)
plot(universal.treeshape)

## Visual representation of two trees at the same time
data(carnivora.treeshape)
plot(carnivora.treeshape, cutreeshape(carnivora.treeshape, 8, "bottom"))
}

\keyword{hplot}
