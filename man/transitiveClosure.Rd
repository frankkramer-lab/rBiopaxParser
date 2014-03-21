\name{transitiveClosure}
\alias{transitiveClosure}
\title{This function generates the transitive closure of the supplied graph.}
\usage{
transitiveClosure(mygraph)
}
\arguments{
  \item{mygraph}{graphNEL}
}
\value{
Returns the transitive closure of the supplied graph.
}
\description{
This function generates the transitive closure of the
supplied graph. In short: if A->B->C then an edge A->C is
added. Edge weights are conserved if possible (in a
hopefully smart way). This is a simple convenience wrapper
for the RBGL function transitive.closure.
}
\author{
Frank Kramer
}

