\name{readBiopax}
\alias{readBiopax}
\title{This function reads in a Biopax .owl file}
\usage{
  readBiopax(file, verbose = TRUE, generateSummary = TRUE,
    generateDF = TRUE)
}
\arguments{
  \item{file}{string. File name}

  \item{verbose}{logical. Output messages about how parsing
  is going and so on.}

  \item{generateSummary}{logical. Generates and attaches
  the summary for the Biopax file to biopax$summary}

  \item{generateDF}{logical. If this is set to FALSE no
  data.frame is generated, only the file is read in to
  biopax$biopaxxml}
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
  \describe{ \item{biopaxxml}{The XML which was read in
  from file.} \item{summary}{The generated summary
  statistic.} \item{df}{The data.frame representing the
  biopax in R} \item{ns_rdf}{RDF Namespace}
  \item{ns_owl}{OWL Namespace} \item{ns_bp}{Biopax
  Namespace} \item{file}{File name} }
}
\author{
  Frank Kramer
}

