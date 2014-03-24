\name{pathway2RegulatoryGraph}
\alias{pathway2RegulatoryGraph}
\title{This function generates the regulatory graph from the activations/inhibitions of a pathway in a biopax model.}
\usage{
pathway2RegulatoryGraph(biopax, pwid, expandSubpathways = TRUE,
  splitComplexMolecules = TRUE, useIDasNodenames = FALSE, verbose = TRUE)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{pwid}{string}

  \item{expandSubpathways}{logical. If TRUE subpathways are
  expanded into this graph, otherwise only this very
  pathway is used.}

  \item{splitComplexMolecules}{logical. If TRUE every
  complex is split up into its components. This leads to
  splitting a single node with name of the complex into
  several nodes with names of the components, these
  components all have identical edges.}

  \item{useIDasNodenames}{logical. If TRUE nodes of the
  graph are named by their molecule IDs instead of using
  the NAME property. This can help with badly
  annotated/formatted databases.}

  \item{verbose}{logical}
}
\value{
Returns the representing the regulatory graph of the
supplied pathway in a node-edge-list graph.
}
\description{
This functions builds a graph from the pathway components
of the supplied pathway. Only instances of class 'control'
are considered, this leads a functinal graph with all edges
either representing activations or inhibitions. No
transports, no translocation, etc. If desired complexes can
be split up into several nodes, this can sometimes lead to
a more complex and cluttered graph. There can not be
multiple edges between 2 nodes. Whenever duplicated edges
are generated (especially by splitting up complexes) a
warning is thrown.
}
\examples{
# load data
 data(biopaxexample)
 pwid1 = "pid_p_100002_wntpathway"
 pwid2 = "pid_p_100146_hespathway"
 mygraph = pathway2RegulatoryGraph(biopax, pwid1)
 plotRegulatoryGraph(mygraph)
}
\author{
Frank Kramer
}

