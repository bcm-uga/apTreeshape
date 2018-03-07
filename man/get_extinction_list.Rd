\name{get_extinction_list}
\alias{get_extinction_list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal BetaAlphaEta function
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

\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}