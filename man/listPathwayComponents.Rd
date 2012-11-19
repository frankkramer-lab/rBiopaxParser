\name{listPathwayComponents}
\alias{listPathwayComponents}
\title{This function lists all pathway components of a given pathway.}
\usage{
  listPathwayComponents(biopax, id)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string. A pathway ID}
}
\value{
  data.frame
}
\description{
  This function returns a (unique) data.frame listing all
  pathway component IDs, names and classes of the supplied
  pathway.
}
\examples{
# load data
 data(biopax2example)
 listPathwayComponents(biopax, id="pid_p_100002_wntpathway")
}
\author{
  Frank Kramer
}

