\name{hivtree.treeshape}
\alias{hivtree.treeshape}

\title{Phylogenetic Tree of 193 HIV-1 Sequences}

\description{
This data set describes an estimated clock-like phylogeny of 193 HIV-1.
group M sequences sampled in the Democratic Republic of Congo. This data is the conversion of the data from the APE package into the class \code{"treeshape"}.
}

\usage{
data(hivtree.treeshape)
}

\format{
\code{hivtree.treeshape} is an object of class \code{"treeshape"}.
}

\source{
This is a data example from Strimmer and Pybus (2001).
}

\references{
Strimmer, K. and Pybus, O. G. (2001) Exploring the demographic history of DNA sequences using the generalized skyline plot. \emph{Molecular Biology and Evolution}, \bold{18}, 2298 -- 2305.
}

\examples{

data("hivtree.treeshape") 
summary(hivtree.treeshape)
plot(hivtree.treeshape)
}

\keyword{datasets}
