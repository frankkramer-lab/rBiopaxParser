\name{internal_getInstanceProperty_df}
\alias{internal_getInstanceProperty_df}
\title{This function returns all properties of the specified type for an instance.}
\usage{
  internal_getInstanceProperty_df(df, id,
    property = "NAME")
}
\arguments{
  \item{df}{data.frame of biopax instances}

  \item{id}{string}

  \item{property}{string. Attention: All properties in
  Biopax Level 2 are all upper case.}
}
\value{
  Returns a character vector with all properties of the
  selected type for this instance. Returns NULL if no
  property of this type is found.
}
\description{
  This function returns all properties of the specified
  type for an instance. By default this function returns
  the NAME property of an instance.
}
\author{
  fkramer
}

