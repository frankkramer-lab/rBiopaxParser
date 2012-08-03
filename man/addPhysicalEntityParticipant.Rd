\name{addPhysicalEntityParticipant}
\alias{addPhysicalEntityParticipant}
\title{This function adds a new physical entity participant.}
\usage{
  addPhysicalEntityParticipant(biopax,
    referencedPhysicalEntityID, ID = NULL)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{referencedPhysicalEntityID}{string. ID the new
  physicalEntity instance to reference here.}

  \item{ID}{string. ID for the physical entity participant.
  If NULL a new ID is generated with prefix
  "physicalEntityParticipant".}
}
\value{
  Returns the biopax model with the added
  physicalEntityParticipant.
}
\description{
  This function adds a new physical entity participant
  instance, which is a placeholder for physicalEntity class
  instances in interactions. This is a convenience function
  to add physicalEntityParticipant instances, internally
  the function addBiopaxInstance is called.
}
\author{
  fkramer
}

