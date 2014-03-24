\name{combineNodes}
\alias{combineNodes}
\title{This function gracefully combines nodes of a regulatory graph.}
\usage{
combineNodes(nodes, graph, newName)
}
\arguments{
  \item{nodes}{vector of node names specifiying which nodes
  to combine.}

  \item{graph}{graphNEL}

  \item{newName}{string. Name of the newly created node
  that will combine the specified nodes.}
}
\value{
Returns a graph with specified nodes removed.
}
\description{
This gracefully combines nodes from a regulatory graph.
This is basically a wrapper for graph::combineNodes(nodes,
graph, newName, collapseFunction=max). If there are
duplicated edges for the nodes, the maximum edgeweight will
be used for the new connection.
}
\examples{
# load data and retrieve wnt pathway
 data(biopaxexample)
}
\author{
Frank Kramer
}

