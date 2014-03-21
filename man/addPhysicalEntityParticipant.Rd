\name{addPhysicalEntityParticipant}
\alias{addPhysicalEntityParticipant}
\title{This function adds a new physical entity participant.}
\usage{
addPhysicalEntityParticipant(biopax, referencedPhysicalEntityID, id = NULL)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{referencedPhysicalEntityID}{string. ID the new
  physicalEntity instance to reference here.}

  \item{id}{string. ID for the physical entity participant.
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
to add physicalEntityParticipant instances, internally the
function addBiopaxInstance is called.
}
\examples{
biopax = createBiopax(level=2)
biopax = addPhysicalEntity(biopax, class="protein", id="p_id1", NAME="protein1")
biopax = addPhysicalEntityParticipant(biopax, "p_id1", id="PEP_p_id1")
biopax = addPhysicalEntity(biopax, class="protein", id="p_id2", NAME="protein2")
biopax = addPhysicalEntityParticipant(biopax, "p_id2", id="PEP_p_id2")
biopax = addBiochemicalReaction(biopax, LEFT=c("PEP_p_id1"), RIGHT=c("PEP_p_id2"), id="biochem_id1")
biopax$dt
}
\author{
fkramer
}

