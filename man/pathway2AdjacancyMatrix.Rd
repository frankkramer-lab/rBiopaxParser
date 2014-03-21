\name{pathway2AdjacancyMatrix}
\alias{pathway2AdjacancyMatrix}
\title{This function generates an adjacency matrix from the activations/inhibitions of a pathway in a biopax model.}
\usage{
pathway2AdjacancyMatrix(biopax, pwid, expandSubpathways = TRUE,
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
Returns the adjacency matrix representing the regulatory
graph of the supplied pathway.
}
\description{
This function internally first calls
pathway2RegulatoryGraph, then converts the regulatory graph
to an adjacency matrix. See pathway2RegulatoryGraph for
more details.
}
\examples{
# load data
 data(biopax2example)
 pwid1 = "pid_p_100002_wntpathway"
 pwid2 = "pid_p_100146_hespathway"
 pathway2AdjacancyMatrix(biopax, pwid1)
}
\author{
Frank Kramer
}

