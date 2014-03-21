\name{listInteractionComponents}
\alias{listInteractionComponents}
\title{This function lists all components of a given interaction.}
\usage{
listInteractionComponents(biopax, id, splitComplexes = TRUE,
  returnIDonly = FALSE)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string. A complex ID}

  \item{splitComplexes}{logical. If TRUE complexes are
  split up into their components and the added to the
  listing.}

  \item{returnIDonly}{logical. If TRUE only IDs of the
  components are returned. This saves tiem for looking up
  names for every single ID.}
}
\value{
data.frame
}
\description{
This function returns a (unique) data.frame listing IDs,
names and classes of all components of the supplied
interaction.
}
\examples{
# load data
 data(biopax2example)
 listInteractionComponents(biopax, id="ex_i_100036_activator_1")
}
\author{
Frank Kramer
}

