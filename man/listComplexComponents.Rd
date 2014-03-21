\name{listComplexComponents}
\alias{listComplexComponents}
\title{This function lists all components of a given complex.}
\usage{
listComplexComponents(biopax, id, returnIDonly = FALSE)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string. A complex ID}

  \item{returnIDonly}{logical. If TRUE only IDs of the
  components are returned. This saves tiem for looking up
  names for every single ID.}
}
\value{
data.frame
}
\description{
This function returns a (unique) data.frame listing all
component IDs, names and classes of the supplied complex.
}
\examples{
# load data
 data(biopax2example)
 listComplexComponents(biopax, id="ex_m_100650")
}
\author{
Frank Kramer
}

