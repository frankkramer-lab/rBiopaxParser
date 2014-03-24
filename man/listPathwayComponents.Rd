\name{listPathwayComponents}
\alias{listPathwayComponents}
\title{This function lists all pathway components of a given pathway.}
\usage{
listPathwayComponents(biopax, id, includeSubPathways = TRUE,
  returnIDonly = FALSE, biopaxlevel = 3)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string. A pathway ID}

  \item{includeSubPathways}{logical. If TRUE the returned
  list will include subpathways and pathwaysteps as well.}

  \item{returnIDonly}{logical. If TRUE only IDs of the
  components are returned. This saves tiem for looking up
  names for every single ID.}

  \item{biopaxlevel}{integer. Set the biopax level here if
  you supply a data.table directly.}
}
\value{
data.frame
}
\description{
This function returns a (unique) data.frame listing all
component IDs, names and classes of the supplied pathway.
}
\examples{
# load data
 data(biopaxexample)
 listPathwayComponents(biopax, id="pid_p_100002_wntpathway")
}
\author{
Frank Kramer
}

