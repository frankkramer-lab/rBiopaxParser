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
  compatible data.frame) to an existing biopax model.
}
\author{
  Frank Kramer
}

