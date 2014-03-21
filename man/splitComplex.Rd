\name{splitComplex}
\alias{splitComplex}
\title{This functions splits up a complex into its components.}
\usage{
splitComplex(biopax, complexid, recursive = TRUE, returnIDonly = FALSE,
  biopaxlevel = 3)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{complexid}{string ID of an complex}

  \item{recursive}{logical}

  \item{returnIDonly}{logical. If TRUE only IDs of the
  components are returned. This saves tiem for looking up
  names for every single ID.}

  \item{biopaxlevel}{integer. Set the biopax level here if
  you supply a data.table directly.}
}
\value{
Returns a character vector with the names of all
subcomponents.
}
\description{
This function looks up the supplied Complex ID and returns
the names of all its components.
}
\examples{
# load data
 data(biopax2example)
 selectInstances(biopax, id="ex_m_100650")
 listInstances(biopax, id="ex_m_100650")
 listComplexComponents(biopax, id="ex_m_100650")
 splitComplex(biopax, complexid="ex_m_100650")
}
\author{
Frank Kramer
}

