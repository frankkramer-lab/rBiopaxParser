\name{listPathways}
\alias{listPathways}
\title{This function returns a list of all pathway ids.}
\usage{
listPathways(biopax, biopaxlevel = 3)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{biopaxlevel}{integer. Set the biopax level here if
  you supply a data.table directly.}
}
\value{
Returns a character vector containing the names of all
pathways.
}
\description{
This function returns a vector of all pathway ids.
}
\examples{
# load data
 data(biopax2example)
 listPathways(biopax)
}
\author{
Frank Kramer
}

