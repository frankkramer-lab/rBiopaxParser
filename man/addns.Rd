\name{addns}
\alias{addns}
\title{Add a namespace tag to the supplied classname string}
\usage{
addns(classname, namespace = "bp")
}
\arguments{
  \item{classname}{A string containing a classname}

  \item{namespace}{A string containing a namespace}
}
\value{
If the classname is not preceeded by a namespace yet, the
supplied namespace is pasted in front of it and returned.
}
\description{
This function takes the input classname, checks if it
already has a namespace, and if not pastes the namespace
tag with a dividing ":" in front of it.
}
\author{
Frank Kramer
}

