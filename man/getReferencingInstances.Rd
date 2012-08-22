\name{getReferencingInstances}
\alias{getReferencingInstances}
\title{This function returns a data.frame containing all instances referencing the supplied instance.}
\usage{
  getReferencingInstances(biopax, id, recursive = TRUE,
    onlyFollowProperties = c())
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string}

  \item{recursive}{logical}

  \item{onlyFollowProperties}{character vector}
}
\value{
  Returns a data.frame containing all referencing
  instances.
}
\description{
  This function takes an id and a biopax model as input. A
  data.frame containing all instances referencing the
  supplied instance is returned. If recursive == TRUE this
  function recurses through all referencing IDs of the
  referencing instances and so on. "onlyFollowProperties"
  limits the recursivness to only certain properties, for
  example follow only complexes or physicalEntities.
}
\author{
  Frank Kramer
}

