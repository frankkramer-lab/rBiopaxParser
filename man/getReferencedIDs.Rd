\name{getReferencedIDs}
\alias{getReferencedIDs}
\title{This function returns a vector of ids of all instances referenced by the specified instance.}
\usage{
getReferencedIDs(biopax, id, recursive = TRUE, onlyFollowProperties = c())
}
\arguments{
  \item{biopax}{A biopax model OR a compatible data.table}

  \item{id}{string. ID of the instance}

  \item{recursive}{logical}

  \item{onlyFollowProperties}{character vector}
}
\value{
Returns a character vector of IDs referenced by the
supplied id in the supplied biopax model.
}
\description{
This function takes an id and a biopax model as input. The
id of every instance that is referenced is returned. If
recursive == TRUE this function recurses through all
referenced IDs of the referenced instances and so on.
"onlyFollowProperties" limits the recursivness to only
certain properties, for example follow only complexes or
physicalEntities.
}
\examples{
# load data
 data(biopaxexample)
 listComplexComponents(biopax, id="ex_m_100650")
 getReferencedIDs(biopax, id="ex_m_100650", recursive=FALSE)
 getReferencedIDs(biopax, id="ex_m_100650", recursive=TRUE)
}
\author{
Frank Kramer
}

