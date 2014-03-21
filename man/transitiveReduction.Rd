\name{transitiveReduction}
\alias{transitiveReduction}
\title{This function generates the transitive reduction of the supplied graph.}
\usage{
transitiveReduction(mygraph)
}
\arguments{
  \item{mygraph}{graphNEL}
}
\value{
Returns the transitive reduction of the supplied graph.
}
\description{
This function generates the transitive reduction of the
supplied graph. In short: if A->B->C AND A->C then edge
A->C is removed. This is a simple convenience wrapper for
the NEM function transitive.reduction. Be aware of
implications on the edge weights!
}
\author{
Frank Kramer
}

