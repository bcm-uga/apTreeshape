\name{get_tree_beta}
\alias{get_tree_beta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Beta parameter as a function of the proportion of remaining tips
}
\description{
Computes the maximum likelihood estimate of the parameter beta in trees simulated by the model, as a function of the proportion of conserved species
}
\usage{
get_tree_beta(epsilon, beta, alpha, N, sampl.frac, ntree, equal.ab = TRUE, eta = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epsilon}{
Minimum size of unsampled splits (see appendix 1)
}
  \item{beta}{
Imbalance index
}
  \item{alpha}{
Clade age-richness index
}
  \item{N}{
Initial tip number
}
  \item{sampl.frac}{
Vector of tips fractions for which we want to compute the beta statistic
}
  \item{ntree}{
Number of simulated trees
}
  \item{equal.ab}{
If set to TRUE, all species have the same probability to go extinct first (default to TRUE)
}
  \item{eta}{
Clade abundance-richness index (if equal.ab == FALSE)
}
}

\value{
A table of size length(sample.fract)*ntree. The element in ligne i and column j is the Maximum Likelihood Estimate for the beta statistic of the j^th tree in wich a fraction sample.frac[i] has been sampled.  
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}
\examples{
# With the field of bullets hypothesis

set.seed(813)
sampl.frac=seq(0.1,1,0.05)
Beta=get_tree_beta(epsilon=0.001, beta=-1, alpha=-1, N=50, sampl.frac=sampl.frac, ntree=10)

Beta_quantiles=sapply(1:nrow(Beta),function(x){quantile(Beta[x,],c(0.05,0.5,0.95))})

plot(1, type="n", xlab="Fraction of extinct species, p", ylab="Beta statistic", 
      ylim=c(-2,10), xlim=c(0,1))
polygon(c(1-sampl.frac, rev(1-sampl.frac)), c(Beta_quantiles[1,(1:length(sampl.frac))],
        rev(Beta_quantiles[3,(1:length(sampl.frac))])), border=NA, col=grey(0.7))
points(1-sampl.frac, Beta_quantiles[2,(1:length(sampl.frac))],t="l")


# With nonrandom extinctions

set.seed(813)
sampl.frac=seq(0.05,1,0.05)
Beta=get_tree_beta(epsilon=0.001, beta=5, alpha=2, eta=2, N=50, 
                    sampl.frac=sampl.frac, ntree=10)

Beta_quantiles=sapply(1:nrow(Beta),function(x){quantile(Beta[x,],c(0.05,0.5,0.95))})

plot(1, type="n", xlab="Fraction of extinct species, p", 
      ylab="Beta statistic", ylim=c(-2,10), xlim=c(0,1))
polygon(c(1-sampl.frac, rev(1-sampl.frac)), 
      c(Beta_quantiles[1,(1:length(sampl.frac))], 
      rev(Beta_quantiles[3,(1:length(sampl.frac))])), border=NA, col=grey(0.7))
points(1-sampl.frac, Beta_quantiles[2,(1:length(sampl.frac))],t="l")


}

