\name{getInstanceProperty}
\alias{getInstanceProperty}
\title{This function returns all properties of the specified type for an instance.}
\usage{
getInstanceProperty(biopax, id, property = "NAME", includeAllNames = TRUE)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string}

  \item{property}{string.}

  \item{includeAllNames}{logical. Biopax Level 3 brought 2
  new name properties: displayName and standardName. Per
  default this return all names of an instance. Disable if
  you only want the NAME property.}
}
\value{
Returns a character vector with all properties of the
selected type for this instance. Returns NULL if no
property data is found.
}
\description{
This function returns all properties of the specified type
for an instance. By default this function returns the NAME
property of an instance.
}
\examples{
# load data
 data(biopax2example)
 getInstanceProperty(biopax, id="ex_m_100650", property="NAME")
 getInstanceProperty(biopax, id="ex_m_100650", property="ORGANISM")
 getInstanceProperty(biopax, id="ex_m_100650", property="COMPONENTS")
}
\author{
fkramer
}

