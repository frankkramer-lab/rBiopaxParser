\name{getXrefAnnotations}
\alias{getXrefAnnotations}
\title{This function returns the annotations of the supplied instances.}
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
\examples{
# load data
 data(biopax2example)
# no annotations for exactly the complex
getXrefAnnotations(biopax, id="ex_m_100650")
# split up the complex and get annotations for all the molecules involved
getXrefAnnotations(biopax, id="ex_m_100650", splitComplexes=TRUE)
}
\author{
  fkramer
}

