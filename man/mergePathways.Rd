\name{mergePathways}
\alias{mergePathways}
\title{This function merges two given pathways}
\usage{
  mergePathways(biopax, pwid1, pwid2, NAME, id = NULL,
    ORGANISM = "", COMMENT = NULL)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{pwid1}{string. ID of first pathway to merge}

  \item{pwid2}{string. ID of second pathway to merge}

  \item{NAME}{string. Name of the new merged pathway}

  \item{id}{string. ID for the pathway. If NULL a new ID is
  generated with prefix "pathway".}

  \item{ORGANISM}{string. Organism property of the pathway.
  By default uses the same organism as the first supplied
  pathway. If NULL no organism property is set.}

  \item{COMMENT}{string. An optional comment}
}
\value{
  A biopax model with the merged pathway added.
}
\description{
  This function merges two given pathways and appends it to
  the supplied biopax model. The user has to specify a new
  name for the pathways and can supply ID, ORGANISM and
  COMMENT properties for the new pathway. If no ID is
  supplied, a new unique ID is generated. If no organism
  property is supplied the organism property of the first
  pathway is re-used. If ORGANISM is NULL the property is
  not set. Optionally a comment can be added to the
  pathway.
}
\author{
  fkramer
}

