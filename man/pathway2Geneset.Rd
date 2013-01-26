\name{pathway2Geneset}
\alias{pathway2Geneset}
\title{This function generates the gene set of a pathway.}
\usage{
  pathway2Geneset(biopax, pwid)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{pwid}{string}
}
\value{
  Returns the gene set of the supplied pathway. Returns
  NULL if the pathway has no components.
}
\description{
  This function generates a gene set of all
  physicalEntity's of a pathway. First all interactions of
  the pathway are retrieved and all components of these
  interactions are then listed.
}
\examples{
# load data
 data(biopax2example)
 pwid1 = "pid_p_100002_wntpathway"
 pathway2Geneset(biopax, pwid=pwid1)
}
\author{
  Frank Kramer
}

