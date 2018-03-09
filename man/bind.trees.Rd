\name{bind.trees}
\alias{bind.trees}
\title{
Binds two tree together
}
\description{
Binds two tree together
}
\usage{
bind.trees(a1, a2, root1, root2, ab = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a1}{
A phylo object
}
  \item{a2}{
A phylo object
}
  \item{root1}{
The length of the branch between the crown of the resulting tree and the crown of a1
}
  \item{root2}{
The length of the branch between the crown of the resulting tree and the crown of a2
}
  \item{ab}{
Boolean indicating wether phylogenies come with tip abundance data (default to FALSE)
}
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}