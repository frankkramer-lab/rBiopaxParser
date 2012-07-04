\name{addPropertiesToBiopaxInstance}
\alias{addPropertiesToBiopaxInstance}
\title{This function adds new properties to an existing biopax instance.}
\usage{
  addPropertiesToBiopaxInstance(biopax, instanceid,
    properties)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{instanceid}{string. ID of the instance}

  \item{properties}{named list of properties.}
}
\value{
  Returns the supplied biopax model with new properties
  added to this instance.
}
\description{
  This function adds new properties to an existing biopax
  instance.
}
\author{
  Frank Kramer
}

