\name{addPathway}
\alias{addPathway}
\title{This function adds a new pathway to the biopax model.}
\usage{
  addPathway(biopax, NAME, PATHWAY_COMPONENTS = c(),
    id = NULL, ORGANISM = NULL, COMMENT = NULL)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{NAME}{string. Name of the pathway}

  \item{PATHWAY_COMPONENTS}{character vector. IDs of the
  pathway components. This must be IDs of instances of type
  interaction/pathway/pathwayStep (or their subclasses).}

  \item{id}{string. ID for the pathway. If NULL a new ID is
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
\examples{
biopax = createBiopax(level=2)
biopax = addPhysicalEntity(biopax, class="protein", id="p_id1", NAME="protein1")
biopax = addPhysicalEntityParticipant(biopax, "p_id1", id="PEP_p_id1")
biopax = addPhysicalEntity(biopax, class="protein", id="p_id2", NAME="protein2")
biopax = addPhysicalEntityParticipant(biopax, "p_id2", id="PEP_p_id2")
biopax = addBiochemicalReaction(biopax, LEFT=c("PEP_p_id1"), RIGHT=c("PEP_p_id2"), id="biochem_id_1")
biopax = addPhysicalEntity(biopax, class="protein", id="p_id3", NAME="controllerProtein1")
biopax = addPhysicalEntityParticipant(biopax, "p_id3", id="PEP_p_id3")
biopax = addControl(biopax, CONTROL_TYPE="ACTIVATION", CONTROLLER="PEP_p_id3", CONTROLLED="biochem_id_1", id="c_id1")
biopax = addPathway(biopax, NAME="mypathway1", PATHWAY_COMPONENTS=c("c_id1"), id="pw_id1")
biopax$df
}
\author{
  fkramer
}

