\name{isOfClass}
\alias{isOfClass}
\title{Checks if instances in the biopax data.table are of the given class}
\usage{
isOfClass(df, class, considerInheritance = FALSE, biopaxlevel = 2)
}
\arguments{
  \item{df}{A data.frame with biopax instances}

  \item{class}{A string containing the class name to check
  for}

  \item{considerInheritance}{Logical value indicating
  wether to consider inheritance or not}

  \item{biopaxlevel}{Numeric. Specifies the Biopax Level to
  use.}
}
\value{
Returns TRUE for every row in the data.frame which is of
the supplied class
}
\description{
This function checks if instances in the supplied biopax
data.table are of a given class. If considerInheritance is
set to TRUE it also checks if instances are of a given
class or any of its inherited classes.
}
\author{
Frank Kramer
}

