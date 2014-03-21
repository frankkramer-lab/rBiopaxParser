\name{addPhysicalEntity}
\alias{addPhysicalEntity}
\title{This function adds a new physical entity.}
\usage{
addPhysicalEntity(biopax, class = c("dna", "rna", "protein", "smallMolecule",
  "complex")[1], NAME, id = NULL, ORGANISM = NULL, COMMENT = NULL)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{class}{string. Class of the physical entity to add,
  choose from
  c("dna","rna","protein","smallMolecule","complex").}

  \item{NAME}{string. Name of the new physical entity}

  \item{id}{string. ID for the physical entity. If NULL a
  new ID is generated with prefix "physicalEntity".}

  \item{ORGANISM}{string. Organism property of the
  molecule. optional.}

  \item{COMMENT}{string. An optional comment}
}
\value{
Returns the biopax model with the added physical entity.
}
\description{
This function adds a new physical entity of chosen class to
the biopax model. This is a convenience function to add
physical entities, internally the function
addBiopaxInstance is called with properties NAME and
ORGANISM set.
}
\examples{
biopax = createBiopax(level=2)
biopax = addBiopaxInstance(biopax, class="protein", id="id1", properties=list(NAME="protein1",COMMENT="this is my first protein!"))
biopax$dt
biopax = addPhysicalEntity(biopax, class="protein", id="id2", NAME="protein2", COMMENT="This is a protein added using the convenience function addPhysicalEntitiy")
biopax$dt
}
\author{
fkramer
}

