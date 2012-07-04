\name{getInstanceProperties}
\alias{getInstanceProperties}
\title{This function returns all properties of the specified type for an instance.}
\usage{
  getInstanceProperties(biopax, instanceid,
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
  selected type for this instance.
}
\description{
  This function returns all properties of the specified
  type for an instance.
}
\author{
  fkramer
}

