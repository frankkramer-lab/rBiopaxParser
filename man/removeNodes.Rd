\name{removeNodes}
\alias{removeNodes}
\title{This function gracefully removes nodes from a regulatory graph.}
\usage{
  removeNodes(graph, nodes)
}
\arguments{
  \item{graph}{graphNEL}

  \item{nodes}{vector of node names specifiying which nodes
  to remove.}
}
\value{
  Returns a graph with specified nodes removed.
}
\description{
  This function gracefully removes nodes from a regulatory
  graph. If the node to be removed has both parent and
  child nodes, these are connected directly. The weight of
  the new direct edge is the product of multiplying the
  incomming and outgoing edge weights of the original node.
}
\examples{
# load data and retrieve wnt pathway
 data(biopax2example)
 pwid1 = "pid_p_100002_wntpathway"
 mygraph1 = pathway2RegulatoryGraph(biopax, pwid1)
 mygraph1 = layoutRegulatoryGraph(mygraph1)
 # retrieve all nodes
 nodes = nodes(mygraph1)
 # random expression data for your nodes
 values = rnorm(length(nodes), mean=6, sd=2)
 # color nodes of the graph
 mygraph1 = colorGraphNodes(mygraph1, nodes, values, colors="greenred")
 # plot the now colored graph
 plotRegulatoryGraph(mygraph1, layoutGraph=FALSE)
}
\author{
  Frank Kramer
}

