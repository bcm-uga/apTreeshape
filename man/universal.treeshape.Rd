\name{universal.treeshape}
\alias{universal.treeshape}

\title{ Universal phylogenetic tree of life }

\description{
This data set describes the Universal Tree of Life described by CR Woese, in the \code{treeshape} class.
"The universal phylogenetic tree not only spans all extant life, but its root and earliest branchings represent stages in the evolutionary process before modern cell types had come into being."
}

\usage{
data(universal.treeshape)
}

\format{
\code{universal.treeshape} is an object of class \code{"treeshape"}.
}

\references{
Woese, C. R. (July, 2000) \emph{Interpreting the universal phylogenetic tree}, \bold{97}, 392 -- 8396 (Proc. Natl. Acad. Sci. USA).
}

\examples{

## Example tree in "treeshape" format
data("universal.treeshape") 

## Summary of the tree
summary(universal.treeshape)

## Visual representation of the tree
plot(universal.treeshape)
}

\keyword{datasets}
