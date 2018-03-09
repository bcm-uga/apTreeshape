\name{get_extinction_list}
\alias{get_extinction_list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gives the tips of the phylogeny in their extinction order
}
\description{
Orders the tips which will sequentially get extinct/be sampled, based on relative species abundance as explained in the main text of the article.
}
\usage{
get_extinction_list(N, tree, equal.ab = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Tip nuber
}
  \item{tree}{
A phylo object
}
  \item{equal.ab}{
If set to TRUE (the default), all species have the same probability to go extinct first}
}
\references{
Maliet O., Gascuel F., Lambert A. (2018) \emph{Ranked tree shapes, non-random 
extinctions and the loss of phylogenetic diversity}, bioRxiv 224295, 
doi: https://doi.org/10.1101/224295 
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}