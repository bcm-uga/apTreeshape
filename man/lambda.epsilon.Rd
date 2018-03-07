\name{lambda.epsilon}
\alias{lambda.epsilon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal BetaAlphaEta function
}
\description{
Computes the lambda_epsilon from appendix 1
}
\usage{
lambda.epsilon(epsilon, beta, n, upper = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epsilon}{
Minimum size of unsampled splits (see appendix 1)
}
  \item{beta}{
Imbalance index
}
  \item{n}{
Tip number
}
  \item{upper}{
Upper bound from the integral (see lambda_n expression, default to 10000)
}
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}
