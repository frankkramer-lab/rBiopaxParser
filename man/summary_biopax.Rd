\name{summary_biopax}
\alias{summary_biopax}
\title{This function generates a summary statistics for the biopax model.}
\usage{
  summary_biopax(object, verbose = TRUE)
}
\arguments{
  \item{object}{A biopax model}

  \item{verbose}{logical}
}
\value{
  Returns the summary for the supplied biopax model.
}
\description{
  This function generates a summary statistics for the
  biopax model. This function can take a while with really
  big Biopax files like NCIs Pathway Interaction Database
  or Reactome. This function is called internally by
  readBiopax if generateSummary == TRUE, so just check
  biopax$summary if it has already been generated.
}
\author{
  Frank Kramer
}

