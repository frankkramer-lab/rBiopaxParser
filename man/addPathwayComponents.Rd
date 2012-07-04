\name{addPathwayComponents}
\alias{addPathwayComponents}
\title{This function adds pathway components to an existing pathway}
\usage{
  addPathwayComponents(biopax, ID,
    PATHWAY_COMPONENTS = c())
}
\arguments{
  \item{biopax}{A biopax model}

  \item{ID}{string. ID for the pathway}

  \item{PATHWAY_COMPONENTS}{character vector. IDs of the
  pathway components. This must be IDs of instances of type
  interaction/pathway/pathwayStep (or their subclasses).}
}
\value{
  Returns the biopax model with the pathway components
  added to the pathway
}
\description{
  This function adds pathway components to an existing
  pathway. Property PATHWAY-COMPONENTS are references to
  IDs of interaction/pathways/pathwaySteps (or subclasses
  of those)
}
\author{
  fkramer
}

