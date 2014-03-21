\name{internal_generateXMLfromBiopax}
\alias{internal_generateXMLfromBiopax}
\title{This function generates the xmlTree from the supplied biopax model.}
\usage{
internal_generateXMLfromBiopax(biopax, namespaces = namespaces,
  verbose = TRUE)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{namespaces}{A list of namespaces to use for the
  generated XML/RDF file}

  \item{verbose}{logical}
}
\value{
Returns the xmlTree generated from the supplied biopax
model.
}
\description{
This function is used internally by writeBiopax. It can
also be called directly with a fitting dataframe in
list(df=data.frame()), but this will probably break things.
}
\author{
Frank Kramer
}

