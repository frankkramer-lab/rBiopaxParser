\name{getReferencedIDs}
\alias{getReferencedIDs}
\title{This function returns a vector of ids of all referenced instances of the supplied instance.}
\usage{
  getReferencedIDs(biopax, id, recursive = TRUE,
    onlyFollowProperties = c())
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string. ID of the instance}

  \item{recursive}{logical}

  \item{onlyFollowProperties}{character vector}
}
\value{
  Returns a character vector of IDs referenced by the
  supplied id in the supplied biopax model.
}
\description{
  This function takes an id and a biopax model as input.
  The id of every instance that is referenced is returned.
  If recursive == TRUE this function recurses through all
  referenced IDs of the referenced instances and so on.
  "onlyFollowProperties" limits the recursivness to only
  certain properties, for example follow only complexes or
  physicalEntities. (not yet implemented)
}
\author{
  Frank Kramer
}

