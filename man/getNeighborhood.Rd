\name{getNeighborhood}
\alias{getNeighborhood}
\title{This function returns the neighborhood of a physicalEntity}
\usage{
  getNeighborhood(biopax, id, depth = 1,
    onlyInPathways = c())
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string. ID of a physicalEntity (dna, rna,
  protein, complex, smallMolecule)}

  \item{depth}{integer. Search depth, this specifies how
  far out from the specified molecule the neighborhood
  should be streched.}

  \item{onlyInPathways}{character vector of pathway IDs.
  Search only in these pathways for neighbors.}
}
\value{
  Returns ids of interactions within 'depth' number of
  steps of the specified physicalEntity
}
\description{
  This function searches the supplied biopax for
  interactions that are connected to the molecule or within
  'depth' number of steps from it.
}
\author{
  fkramer
}

