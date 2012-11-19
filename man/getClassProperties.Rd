\name{getClassProperties}
\alias{getClassProperties}
\title{This function returns the properties of the supplied biopax class.}
\usage{
  getClassProperties(classname)
}
\arguments{
  \item{classname}{A string containing a class name}
}
\value{
  Returns character vector containing the superclasses of
  the supplied class
}
\description{
  This function returns the superclasses of the supplied
  biopax class. It always considers inhertance. Every class
  inhertis the properties of its super classes. A table
  listing all available properties and their cardinalities
  (for Biopax Level 2).
}
\examples{
getClassProperties("control")
}
\author{
  Frank Kramer
}

