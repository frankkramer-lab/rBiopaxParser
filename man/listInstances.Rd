\name{listInstances}
\alias{listInstances}
\title{Lists all instances that conform to the selection criteria.}
\usage{
listInstances(biopax, id = NULL, class = NULL, name = NULL,
  includeSubClasses = FALSE, returnIDonly = FALSE, biopaxlevel = 3)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string. ID of the instances to select}

  \item{class}{string. Class of the instances to select}

  \item{name}{string. Name of the instances to select}

  \item{includeSubClasses}{logical. If includeSubClasses is
  set to TRUE the class criteria is broadened to include
  all classes that inherit from the given class}

  \item{returnIDonly}{logical. If TRUE only IDs of the
  components are returned. This saves time for looking up
  names for every single ID.}

  \item{biopaxlevel}{integer. Set the biopax level here if
  you supply a data.table directly.}
}
\value{
Returns a data.frame containing all instances conforming to
the given selection criteria. If returnIDonly=TRUE, only
the selector for the internal data.table otherwise.
}
\description{
Lists all instances that conform to the selection criteria.
In contrast to selectInstances this function returns an
easier to read list. This function returns an ordered
data.table of class, id and name of the instances.
Selection criteria are wether instances belong to a certain
class or have the specified id or name. Setting a criteria
to NULL ignores this criteria. If includeSubClasses is set
to TRUE the class criteria is broadened to include all
classes that inherit from the given class, e.g. if
class="control" and includeSubClasses=TRUE the function
will select catalyses and modulations too, since they are a
subclass of class control.
}
\examples{
# load data
 data(biopaxexample)
 # list all instances of class "protein"
 listInstances(biopax, class="protein")
 # list all instances of class "pathway"
 listInstances(biopax, class="pathway")
 # list all interaction including all subclasses of interactions
 listInstances(biopax, class="interaction", includeSubClasses=TRUE)
}
\author{
Frank Kramer
}

