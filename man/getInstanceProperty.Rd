\name{getInstanceProperty}
\alias{getInstanceProperty}
\title{This function returns all properties of the specified type for an instance.}
\usage{
  getInstanceProperty(biopax, instanceid,
    property = "NAME")
}
\arguments{
  \item{biopax}{A biopax model}

  \item{instanceid}{string}

  \item{property}{string. Attention: All properties in
  Biopax Level 2 are all upper case.}
}
\value{
  Returns a character vector with all properties of the
  selected type for this instance. Returns NULL if no
  property of this type is found.
}
\description{
  This function returns all properties of the specified
  type for an instance.
}
\author{
  fkramer
}

