\name{getClassProperties}
\alias{getClassProperties}
\title{This function returns the properties of the supplied biopax class.}
\usage{
getClassProperties(classname, biopaxlevel = 3)
}
\arguments{
  \item{classname}{A string containing a class name}

  \item{biopaxlevel}{Numeric. Specifies the Biopax Level to
  use.}
}
\value{
Returns a data.frame containing the properties and
cardinalities of the supplied class
}
\description{
This function returns the properties of the supplied biopax
class. It always considers inhertance. Every class inhertis
the properties of its super classes. A table listing all
available properties and their cardinalities (for Biopax
Level 2).
}
\examples{
getClassProperties("control")
}
\author{
Frank Kramer
}

