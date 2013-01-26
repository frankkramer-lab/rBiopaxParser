\name{getSuperClasses}
\alias{getSuperClasses}
\title{This function returns the superclasses of the supplied biopax class.}
\usage{
  getSuperClasses(classname, biopaxlevel = 2)
}
\arguments{
  \item{classname}{A string containing a class name}

  \item{biopaxlevel}{Numeric. Specifies the Biopax Level to
  use.}
}
\value{
  Returns character vector containing the superclasses of
  the supplied class
}
\description{
  This function returns the superclasses of the supplied
  biopax class.
}
\examples{
getSuperClasses("control")
}
\author{
  Frank Kramer
}

