\name{getReferencingIDs}
\alias{getReferencingIDs}
\title{This function returns a vector of ids of all instances that reference the supplied id.}
\usage{
  getReferencingIDs(biopax, id, recursive = TRUE,
    onlyFollowProperties = c())
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string. ID of the instance}

  \item{recursive}{logical}

  \item{onlyFollowProperties}{character vector}
}
\value{
  Returns a character vector of IDs referencing the
  supplied id in the supplied biopax model.
}
\description{
  This function takes an id and a biopax model as input.
  The id of every instance that references the supplied id
  is returned. If recursive == TRUE this function recurses
  through all referencing IDs of the referencing instances
  and so on. "onlyFollowProperties" limits the recursivness
  to only certain properties, for example follow only
  complexes or physicalEntities.
}
\examples{
# load data
 data(biopax2example)
 listComplexComponents(biopax, id="ex_m_100650")
 getReferencingIDs(biopax, id="ex_m_100650", recursive=FALSE)
 getReferencingIDs(biopax, id="ex_m_100650", recursive=TRUE)
}
\author{
  Frank Kramer
}

