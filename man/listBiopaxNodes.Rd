\name{listBiopaxNodes}
\alias{listBiopaxNodes}
\title{This function returns a summary statistic for the biopax model.}
\usage{
  listBiopaxNodes(biopax)
}
\arguments{
  \item{biopax}{A biopax model}
}
\value{
  This function returns a summary statistic for the biopax
  model.
}
\description{
  This function returns a summary statistic for the biopax
  model. It checks if the summary has already been
  generated with summary.biopax and generates it if
  necessary. Have a look at the biopax$summary entry for
  more statistics. This function can take a while with
  really big Biopax files like NCIs Pathway Interaction
  Database or Reactome.
}
\author{
  Frank Kramer
}

