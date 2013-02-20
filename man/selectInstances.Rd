\name{selectInstances}
\alias{selectInstances}
\title{Returns all instances that conform to the selection criteria.}
\usage{
  selectInstances(biopax, id = NULL, class = NULL,
    property = NULL, name = NULL, returnValues = TRUE,
    includeSubClasses = FALSE,
    includeReferencedInstances = FALSE)
}
\arguments{
  \item{biopax}{A biopax model}

  \item{id}{string. ID of the instances to select}

  \item{class}{string. Class of the instances to select}

  \item{property}{string. Return only this property of the
  instances}

  \item{name}{string. Name of the instances to select}

  \item{returnValues}{logical. If returnValues is set to
  FALSE only the selector (a logical vector with length of
  the internal data.frame) is returned, otherwise the
  selected data is returned}

  \item{includeSubClasses}{logical. If includeSubClasses is
  set to TRUE the class criteria is broadened to include
  all classes that inherit from the given class}

  \item{includeReferencedInstances}{logical. If
  includeReferencedInstances is set to TRUE all instances
  that are being referenced by the selected instances are
  being selected too}
}
\value{
  Returns a data.frame containing all instances conforming
  to the given selection criteria if returnValues=TRUE,
  only the selector for the internal data.frame otherwise.
}
\description{
  Returns all instances that conform to the selection
  criteria. This function returns a subset of the internal
  data.frame of the biopax object. Selection criteria are
  wether instances belong to a certain class or have the
  specified id, property or name. Setting a criteria to
  NULL ignores this criteria. If returnValues is set to
  FALSE only the selector (a logical vector with length of
  the internal data.frame) is returned, otherwise the
  selected data is returned. If includeSubClasses is set to
  TRUE the class criteria is broadened to include all
  classes that inherit from the given class, e.g. if
  class="control" and includeSubClasses=TRUE the function
  will select catalyses and modulations too, since they are
  a subclass of class control. If
  includeReferencedInstances is set to TRUE all instances
  that are being referenced by the selected instances are
  being selected too. The parameter works recursively, this
  means for example that a selected pathway and all it's
  interactions, complexes, molecules and annotations are
  returned if this parameter is set to true. This parameter
  is especially helpful if you want to migrate or merge
  knowledge from different data bases.
}
\examples{
# load data
 data(biopax2example)
 # select the subset of the internal data.frame that belongs to class "protein"
 selectInstances(biopax, class="protein")
 # select the subset of the internal data.frame that belongs to class "interaction"
 selectInstances(biopax, class="interaction")
 # select the subset of the internal data.frame that belongs to class "interaction" or any of its sub classes, like control, catalysis etc.
 selectInstances(biopax, class="interaction", includeSubClasses=TRUE)
 # select the subset of the internal data.frame that belongs to class "pathway" AND is a "NAME" property
 selectInstances(biopax, class="pathway", property="NAME")
}
\author{
  Frank Kramer
}

