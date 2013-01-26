\name{getSubClasses}
\alias{getSubClasses}
\title{This function returns the subclasses of the supplied biopax class.}
\usage{
  getSubClasses(classname, biopaxlevel = 2)
}
\arguments{
  \item{classname}{A string containing a class name}

  \item{biopaxlevel}{Numeric. Specifies the Biopax Level to
  use.}
}
\value{
  Returns character vector containing the subclasses of the
  supplied class
}
\description{
  This function returns the subclasses of the supplied
  biopax class.
}
\examples{
getSubClasses("control")
}
\author{
  Frank Kramer
}

