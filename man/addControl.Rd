\name{addControl}
\alias{addControl}
\title{This function adds a new control to the biopax model.}
\usage{
  addControl(biopax,
    CONTROL_TYPE = c("ACTIVATION", "INHIBITION"),
    CONTROLLER = "", CONTROLLED = c(), ID = NULL)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{CONTROL_TYPE}{string. Specifies wether this is an
  activating or inhibiting control.}

  \item{CONTROLLER}{string. ID of the
  physicalEntityParticipant instance that is the controller
  of this interaction.}

  \item{CONTROLLED}{vector of strings. IDs of the
  interaction and/or pathway instances that are being
  controlled.}

  \item{ID}{string. ID for the control. If NULL a new ID is
  generated with prefix "control".}
}
\value{
  Returns the biopax model with the added pathway.
}
\description{
  This function adds a new interaction of class control to
  the biopax model. This is a convenience function to add
  controls, internally the function addBiopaxInstance is
  called with properties CONTROL-TYPE, CONTROLLER and
  CONTROLLED set.
}
\author{
  fkramer
}

