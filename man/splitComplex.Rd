\name{splitComplex}
\alias{splitComplex}
\title{This functions splits up a complex into its components.}
\usage{
  splitComplex(biopax, complexid, recursive = TRUE)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{complexid}{string ID of an complex}

  \item{recursive}{logical}
}
\value{
  Returns a character vector with the names of all
  subcomponents.
}
\description{
  This function looks up the supplied Complex ID and
  returns the names of all its components.
}
\author{
  Frank Kramer
}

