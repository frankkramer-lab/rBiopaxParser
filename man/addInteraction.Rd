\name{addInteraction}
\alias{addInteraction}
\title{This function adds a new interaction to the biopax model.}
\usage{
  addInteraction(biopax, class = "control",
    CONTROL_TYPE = c("ACTIVATION", "INHIBITION"),
    CONTROLLER = "", CONTROLLED = c(), ID = NULL)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{class}{string. Class of the interaction to add.
  Suggests "control" or its sibilings "catalysis" or
  "modulation".}

  \item{CONTROL_TYPE}{character vector. IDs of the pathway
  components. This must be IDs of instances of type
  interaction/pathway/pathwayStep (or their subclasses).}

  \item{CONTROLLER}{string. Organism property of the
  pathway. optional.}

  \item{CONTROLLED}{string. An optional comment}

  \item{ID}{string. ID for the pathway. If NULL a new ID is
  generated with prefix "pathway".}
}
\value{
  Returns the biopax model with the added pathway.
}
\description{
  This function adds a new interaction of class control to
  the biopax model. This is a convenience function to add
  controls, internall the function addBiopaxInstance is
  called with properties CONTROL_TYPE, CONTROLLER and
  CONTROLLED set.
}
\author{
  fkramer
}

