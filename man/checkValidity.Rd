\name{checkValidity}
\alias{checkValidity}
\title{This function checks the supplied biopax model for validity.}
\usage{
  checkValidity(biopax)
}
\arguments{
  \item{biopax}{A biopax model}
}
\value{
  logical. Returns TRUE is the biopax model is valid Biopax
  Level 2, or FALSE otherwise.
}
\description{
  This function checks the supplied biopax model for
  validity, concerning classes, properties, etc. Not yet
  implemented. Called internally by writeBiopax.
}
\author{
  Frank Kramer
}

