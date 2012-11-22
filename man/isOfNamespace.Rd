\name{isOfNamespace}
\alias{isOfNamespace}
\title{Check if a classname is preceeded by a certain namespace tag like in "namespace:classname"}
\usage{
  isOfNamespace(classname, namespace = "bp")
}
\arguments{
  \item{classname}{A string containing the classname to
  check}

  \item{namespace}{A string giving the namespace to check
  for}
}
\value{
  This function returns TRUE if the supplied classname
  string is preceeded with the supplied namespace string,
  and FALSE if not.
}
\description{
  This function checks if the supplied input string starts
  with a supplied namespace tag
}
\author{
  Frank Kramer
}

