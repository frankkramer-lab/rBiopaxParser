\name{addBiochemicalReaction}
\alias{addBiochemicalReaction}
\title{This function adds a new biochemical reaction to the biopax model.}
\usage{
  addBiochemicalReaction(biopax, LEFT = c(), RIGHT = c(),
    ID = NULL)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{LEFT}{vector of strings. IDs of the
  physicalEntityParticipant instances that are on the left
  side of this reaction.}

  \item{RIGHT}{vector of strings. IDs of the
  physicalEntityParticipant instances that are on the right
  side of this reaction.}

  \item{ID}{string. ID for the control. If NULL a new ID is
  generated with prefix "biochemicalReaction".}
}
\value{
  Returns the biopax model with the added pathway.
}
\description{
  This function adds a new biochemical reaction of class
  biochemicalReaction to the biopax model. This is a
  convenience function, internally the function
  addBiopaxInstance is called with properties LEFT and
  RIGHT set.
}
\author{
  fkramer
}

