\name{internal_propertyListToDF}
\alias{internal_propertyListToDF}
\title{Internal function to build a data.frame from the list of properties for a new instance}
\usage{
  internal_propertyListToDF(class, id, properties,
    namespace_rdf = "rdf")
}
\arguments{
  \item{class}{string. Class name}

  \item{id}{string. ID of the instance}

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
  Frank Kramer
}

