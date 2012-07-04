\name{getPathwayComponents}
\alias{getPathwayComponents}
\title{This function returns a data.frame containing all pathway components.}
\usage{
  getPathwayComponents(biopax, id)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string. Pathway ID}
}
\value{
  Returns a data.frame containing all pathway components.
}
\description{
  This function returns a data.frame containing all
  instances referenced in a pathways PATHWAY-COMPONENTS
  property.
}
\author{
  Frank Kramer
}

