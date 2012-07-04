\name{getBiopaxInstancesList}
\alias{getBiopaxInstancesList}
\title{Returns all list of all instances of a certain class, or all instances.}
\usage{
  getBiopaxInstancesList(biopax, type = "all")
}
\arguments{
  \item{biopax}{A biopax model}

  \item{type}{string. Class of the instances}
}
\value{
  Returns a character vector containing all instances of
  the supplied class.
}
\description{
  Returns a character vector of IDs of all instances of a
  certain class, or all instances if type=="all".
}
\author{
  Frank Kramer
}

