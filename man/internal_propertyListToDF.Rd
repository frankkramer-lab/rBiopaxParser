\name{internal_propertyListToDF}
\alias{internal_propertyListToDF}
\title{Internal function to build a data.frame from the list of properties for a new instance}
\usage{
  internal_propertyListToDF(instancetype, instanceid,
    properties, namespace_rdf = "rdf")
}
\arguments{
  \item{instancetype}{string. Class name}

  \item{instanceid}{string. ID of the instance}

  \item{properties}{named list of properties.}

  \item{namespace_rdf}{string. This defines the rdf
  namespace to use.}
}
\value{
  Returns a data.frame with the new properties for the
  given instance
}
\description{
  Internal function to build a data.frame from the list of
  properties for a new instance
}
\author{
  Frank Kramer #
}

