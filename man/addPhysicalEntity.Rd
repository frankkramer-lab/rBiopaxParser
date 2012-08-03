\name{addPhysicalEntity}
\alias{addPhysicalEntity}
\title{This function adds a new physical entity.}
\usage{
  addPhysicalEntity(biopax,
    class = c("dna", "rna", "protein", "smallMolecule", "complex")[1],
    NAME, ID = NULL, ORGANISM = NULL, COMMENT = NULL)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{class}{string. Class of the physical entity to add,
  choose from
  c("dna","rna","protein","smallMolecule","complex").}

  \item{NAME}{string. Name of the new physical entity}

  \item{ID}{string. ID for the physical entity. If NULL a
  new ID is generated with prefix "physicalEntity".}

  \item{ORGANISM}{string. Organism property of the
  molecule. optional.}

  \item{COMMENT}{string. An optional comment}
}
\value{
  Returns the biopax model with the added physical entity.
}
\description{
  This function adds a new physical entity of chosen class
  to the biopax model. This is a convenience function to
  add physical entities, internally the function
  addBiopaxInstance is called with properties NAME and
  ORGANISM set.
}
\author{
  fkramer
}

