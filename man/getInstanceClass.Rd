\name{getInstanceClass}
\alias{getInstanceClass}
\title{This function returns the class name of the instance.}
\usage{
getInstanceClass(biopax, id)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string}
}
\value{
Returns the class name of the biopax instance.
}
\description{
This function returns the class name of the instance.
}
\examples{
# load data
 data(biopax2example)
 getInstanceClass(biopax, id="ex_m_100650")
}
\author{
fkramer
}

