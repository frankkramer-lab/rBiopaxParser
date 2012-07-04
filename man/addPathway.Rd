\name{addPathway}
\alias{addPathway}
\title{This function adds a new pathway to the biopax model.}
\usage{
  addPathway(biopax, NAME, PATHWAY_COMPONENTS = c(),
    ID = NULL, ORGANISM = NULL, COMMENT = NULL)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{NAME}{string. Name of the pathway}

  \item{PATHWAY_COMPONENTS}{character vector. IDs of the
  pathway components. This must be IDs of instances of type
  interaction/pathway/pathwayStep (or their subclasses).}

  \item{ID}{string. ID for the pathway. If NULL a new ID is
  generated with prefix "pathway".}

  \item{ORGANISM}{string. Organism property of the pathway.
  optional.}

  \item{COMMENT}{string. An optional comment}
}
\value{
  Returns the biopax model with the added pathway.
}
\description{
  This function adds a new pathway + its PATHWAY-COMPONENTS
  (references to interaction/pathways/pathwaySteps)
}
\author{
  fkramer
}

