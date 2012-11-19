\name{readBiopax}
\alias{readBiopax}
\title{This function reads in a Biopax .owl file}
\usage{
  readBiopax(file, verbose = TRUE)
}
\arguments{
  \item{file}{string. File name}

  \item{verbose}{logical. Output messages about how parsing
  is going and so on.}
}
\value{
  A biopax model
}
\description{
  This function reads in a Biopax .owl file and generates
  the internal data.frame format used in this package. This
  function can take a while with really big Biopax files
  like NCIs Pathway Interaction Database or Reactome. In
  almost every case this is your starting point. Returns a
  biopax model, which is a list with named elements:
  \describe{ \item{df}{The data.frame representing the
  biopax in R} \item{ns_rdf}{RDF Namespace}
  \item{ns_owl}{OWL Namespace} \item{ns_bp}{Biopax
  Namespace} \item{file}{File name} }
}
\examples{
\dontrun{biopax = readBiopax(file="biopaxmodel.owl")}
 \dontrun{biopax}
}
\author{
  Frank Kramer
}

