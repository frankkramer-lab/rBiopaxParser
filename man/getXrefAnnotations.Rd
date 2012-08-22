\name{getXrefAnnotations}
\alias{getXrefAnnotations}
\title{This function returns the annotations of the supplied IDs}
\usage{
  getXrefAnnotations(biopax, id, splitComplexes = FALSE,
    followPhysicalEntityParticipants = TRUE)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{vector of strings. IDs of instances to get
  annotations}

  \item{splitComplexes}{logical. If TRUE complexes are
  split up into their components and the annotation of the
  components is added.}

  \item{followPhysicalEntityParticipants}{logical. If TRUE
  physicalEntityParticipants are resolved to their
  corresponding physicalEntities and their annotation is
  added.}
}
\value{
  Returns data.frame with annotations
}
\description{
  This function returns the annotations of the supplied IDs
  in a data.frame.
}
\author{
  fkramer
}

