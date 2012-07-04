\name{getReferencedInstances}
\alias{getReferencedInstances}
\title{This function returns a data.frame containing all referenced instances of the supplied instance.}
\usage{
  getReferencedInstances(biopax, id, recursive = TRUE,
    onlyFollowProperties = c())
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string}

  \item{recursive}{logical}

  \item{onlyFollowProperties}{character vector}
}
\value{
  Returns a data.frame containing all referenced instances.
}
\description{
  This function takes an id and a biopax model as input. A
  data.frame containing all instances referenced by the
  supplied instance is returned. If recursive == TRUE this
  function recurses through all referenced IDs of the
  referenced instances and so on. "onlyFollowProperties"
  limits the recursivness to only certain properties, for
  example follow only complexes or physicalEntities. (not
  yet implemented)
}
\author{
  Frank Kramer
}

