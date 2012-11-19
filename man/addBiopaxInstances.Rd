\name{addBiopaxInstances}
\alias{addBiopaxInstances}
\title{This function adds new instances to an existing biopax model.}
\usage{
  addBiopaxInstances(biopax, newInstancesDF)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{newInstancesDF}{data.frame. Compatible with
  internal Biopax Level 2 implementation.}
}
\value{
  Returns the supplied biopax model with the new instances
  added.
}
\description{
  This function adds new instances (supplied as a
  compatible data.frame) to an existing biopax model via
  rbind. Usually you want to start out at createBiopax and
  addPhysicalEntity and work your way up the ontology
  ladder.
}
\examples{
# load data
 data(biopax2example)
 biopax_temp = createBiopax(level=2)
 biopax_temp = addBiopaxInstance(biopax_temp, class="protein", id="id1", properties=list(NAME="protein1",SYNONYMS="p1"))
 selectInstances(biopax_temp)
 biopax = addBiopaxInstances(biopax, selectInstances(biopax_temp))
}
\author{
  Frank Kramer
}

