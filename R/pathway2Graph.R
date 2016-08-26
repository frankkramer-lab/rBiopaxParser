###############################################################################
#' This function generates a directed graph from all the interactions of a specified pathway in a biopax model. Edges with no direction are indicated by a 0 weight.
#'  
#' @param biopax A biopax model
#' @param pwid string
#' @param expandSubpathways logical. If TRUE subpathways are expanded into this graph, otherwise only this very pathway is used.
#' @param splitComplexMolecules logical. If TRUE every complex is split up into its components. This leads to splitting a single node with name of the complex into several nodes with names of the components, these components all have identical edges. Default value is FALSE
#' @param useIDasNodenames logical. If TRUE nodes of the graph are named by their molecule IDs instead of using the NAME property. This can help with badly annotated/formatted databases.
#' @param verbose logical 
#' @param withDisconnectedParts logical. If TRUE the pathway graph is returned as such, else only the largest connected component is given back
#' @return Returns the a graph object of the specified pathway. Edges with no direction are indicated by a 0 weight.
#' @author Nirupama Benis
#' @export
#' @import data.table
#' @examples
#'  # load data
#' example data

pathway2Graph <- function (biopax, pwid, expandSubpathways = TRUE, splitComplexMolecules = FALSE, useIDasNodenames = TRUE, verbose = FALSE, withDisconnectedParts = TRUE) { 
  
  if (!require(graph)) {
    message(paste("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!", "\n"))
    return(NULL)
  }
  biopaxlevel <- biopax$biopaxlevel
  pwid <- unique(striphash(pwid))
  mygraph <- new(getClassDef("graphNEL", package = "graph"))
  graph::edgemode(mygraph) <- "directed"
  pwComponentList <- listPathwayComponents(biopax, pwid, returnIDonly = TRUE)
  if (length(pwComponentList) == 0) {
    warning("Pathway seems to have no pathway components")
    return(NULL)
  }
  pwComponentList <- selectInstances(biopax, id = pwComponentList, includeReferencedInstances = TRUE, returnCopy = TRUE)
  pwComponentList$property <- tolower(pwComponentList$property)
  setkeyv(pwComponentList, cols <- c("id", "class", "property"))
  pwInteraction <- pwComponentList[tolower(pwComponentList$class) %chin% c("conversion", "biochemicalreaction", "complexassembly", "transport", "degradation", "transportwithbiochemicalreaction", "molecularinteraction", "control", "catalysis", "templatereactionregulation", "modulation", "geneticinteraction", "templatereaction"), ]
  if (length(pwInteraction$id) == 0) {
    warning("warning: pathway2RegulatoryGraph: supplied graph has no pathway components. Returning NULL.")
    return(NULL)
  } else {
    if (verbose) {
      message(paste("Found", length(unique(pwInteraction$id)), "pathway components. Putting them together..."))
    }
  }
  for (i in unique(pwInteraction$id)) {
    instance = pwInteraction[id == i, ]
    nestedInteraction <- FALSE
    type <- NULL
    if (any(isOfClass(instance, c("Conversion", "conversion"), considerInheritance = TRUE, biopaxlevel = biopaxlevel) || grepl("Degradation", instance$class, ignore.case = TRUE) || grepl("degradation", instance$class, ignore.case = TRUE))){ 
      if (biopaxlevel == 2) {
        type <- "left-right"
      } else { 
        type <- as.character(instance[property == "conversiondirection"]$property_value)
      }
      leftParticipantIds <- striphash(as.character(unique(instance[property == "left"]$property_attr_value)))
      rightParticipantIds <- striphash(as.character(unique(instance[property == "right"]$property_attr_value)))
    }
    if (any(isOfClass(instance, c("Control", "control"), considerInheritance = TRUE, biopaxlevel = biopaxlevel) || grepl("TemplateReactionRegulation", instance$class, ignore.case = TRUE) || grepl("templatereactionregulation", instance$class, ignore.case = TRUE))){ 
      if (biopaxlevel == 2) {
        type <- as.character(instance[property == "control-type"]$property_value)
      } else { 
        type <- as.character(instance[property == "controltype"]$property_value)
      }
      leftParticipantIds <- striphash(as.character(unique(instance[property == "controller"]$property_attr_value)))
      rightParticipantIds <- striphash(as.character(unique(instance[property == "controlled"]$property_attr_value)))
    }
    if (length(type) == 0) {
      if (tolower(as.character(instance[1, class])) == "catalysis") {
        type = "ACTIVATION"
      } else { 
        next
      }
    }
    
    leftParticipants <- NA
    for (i2 in leftParticipantIds) { 
      lInstance <- pwComponentList[id == i2, ]
      if (biopax$biopaxlevel == 2) {
        lInstance <- pwComponentList[id == striphash(lInstance[property == "physical-entity"]$property_attr_value), ]
      }
      leftParticipants <- c(leftParticipants, getParticipants(pwComponentList, lInstance, biopaxlevel))
      if (any(isOfClass(lInstance, c("Interaction"), considerInheritance = TRUE, biopaxlevel = biopaxlevel)))
        warning (paste(pwid, i2, ": Left is an Interaction\n"))
    }
    
    rightParticipants <- NA
    for (i2 in rightParticipantIds) { 
      rInstance <- pwComponentList[id == i2, ]
      if (any(isOfClass(rInstance, c("Conversion", "conversion"), considerInheritance = TRUE, biopaxlevel = biopaxlevel)) || 
            any(isOfClass(rInstance, c("TemplateReaction", "templatereaction"), biopaxlevel = biopaxlevel))) {
        leftrights <- striphash(rInstance[property == "left" | property == "right" | property == "product"]$property_attr_value)
        for (i3 in leftrights) {
          leftrightInstance <- pwComponentList[id == i3, ]
          if (biopax$biopaxlevel == 2) {
            leftrightInstance <- pwComponentList[id == striphash(leftrightInstance[property == "physical-entity"]$property_attr_value), ]
          }
          rightParticipants <- c(rightParticipants, getParticipants(pwComponentList, leftrightInstance, biopaxlevel))      
        }
      } else if (any(isOfClass(rInstance, c("Control", "control"), considerInheritance = TRUE, biopaxlevel = biopaxlevel))) {
        nestedInteraction <- TRUE
        rInteractionInstance <- pwComponentList[id == i2, ]
        if (biopax$biopaxlevel == 2) {
          rInteractionInstance <- pwComponentList[id == striphash(rInteractionInstance[property == "physical-entity"]$property_attr_value), ]
        }
        leftNested <- striphash(rInteractionInstance[property == "controller"]$property_attr_value)
        rightParticipants <- c(rightParticipants, striphash(rInteractionInstance[property == "controller"]$property_attr_value))
        rightNestedIds <- striphash(rInteractionInstance[property == "controlled"]$property_attr_value)
        rightNested <- NA
        for (i3 in rightNestedIds){
          rNestedInstance <- pwComponentList[id == i3,]
          if (any(isOfClass(rNestedInstance, c("Conversion", "control"), considerInheritance = TRUE, biopaxlevel = biopaxlevel)) || 
                any(isOfClass(rNestedInstance, c("TemplateReaction", "templatereaction"), biopaxlevel = biopaxlevel))) { 
            leftrightsNested <- striphash(rNestedInstance[property == "left" | property == "right" | property == "product"]$property_attr_value)
            for(i4 in leftrightsNested){ 
              leftrightsNestedInstance <- pwComponentList[id == i4, ]
              if (biopax$biopaxlevel == 2) {
                leftrightsNestedInstance <- pwComponentList[id == striphash(leftrightsNestedInstance[property == "physical-entity"]$property_attr_value), ]
              }
              rightNested <- c(rightNested, getParticipants(pwComponentList, leftrightsNestedInstance, biopaxlevel))
              rightParticipants <- c(rightParticipants, getParticipants(pwComponentList, leftrightsNestedInstance, biopaxlevel))
            }     
          }
        }
      } else {
        if (biopax$biopaxlevel == 2) {
          rInstance <- pwComponentList[id == striphash(rInstance[property == "physical-entity"]$property_attr_value), ]
        }
        rightParticipants <- c(rightParticipants, getParticipants(pwComponentList, rInstance, biopaxlevel))
      }
    }
    
    leftParticipants <- striphash(unique(leftParticipants))
    rightParticipants <- striphash(unique(rightParticipants))
    leftParticipants <- leftParticipants[!is.na(leftParticipants) & !is.null(leftParticipants) & nchar(leftParticipants) > 0]
    rightParticipants <- rightParticipants[!is.na(rightParticipants) & !is.null(rightParticipants) & nchar(rightParticipants) > 0]
    if (nestedInteraction)
      rightNested <- rightNested[!is.na(rightNested)]
    if (length(leftParticipants) == 0 | length(rightParticipants) == 0) {
      next
    }
    if (verbose) {
      message(paste("Adding to graph: ", unique(as.character(instance[, "class"])), "-", pwComponentList$class[pwComponentList$id == i][1], type, ", LeftParticipants: ", paste(leftParticipants, collapse = " "), "RightParticipants: ", paste(rightParticipants, collapse = " ")))
    }
    
    newnodes <- unique(c(leftParticipants, rightParticipants))
    newnodes <- newnodes[!(newnodes %in% graph::nodes(mygraph))]
    if (length(newnodes) > 0) 
      mygraph <- graph::addNode(newnodes, mygraph)
    if (tolower(type) == "inhibition") {
      weight <- -1
    }  else if (tolower(type) == "spontaneous") {
      weight <- 0
    } else {
      weight <- 1
    }
    for (iLefts in leftParticipants) { 
      for (iRights in rightParticipants) { 
        if (iLefts != iRights) {
          if (grepl("spontaneous", tolower(type))) {
            mygraph <- graph::addEdge(from = iLefts, to = iRights, graph = mygraph, weights = weight)
            mygraph <- graph::addEdge(from = iRights, to = iLefts, graph = mygraph, weights = weight)
          }
          if (!(iLefts %in% unlist(graph::edges(mygraph, which = iRights)))) {
            if(nestedInteraction && (iRights %in% rightNested)) {
              mygraph <- graph::addEdge(from = leftNested, to = iRights, graph = mygraph, weights = weight)
              if(verbose) { 
                print("Nested Interaction")
              }
            } else {
              mygraph <- graph::addEdge(from = iLefts, to = iRights, graph = mygraph, weights = weight)
            }
          } else {
            warning(paste("Problem while adding edge (weight = ", weight, ") to graph: Edge already exists between ", iLefts, " and ", iRights, " (weight = ", ifelse(is.na(unlist(graph::edgeWeights(mygraph))[paste(iLefts, iRights, sep = ".")]), "(opposite direction)", unlist(graph::edgeWeights(mygraph))[paste(iLefts, iRights, sep = ".")]), "). Skipping this edge.", sep = ""))
          }
        }
      }
    }
  }
  if (withDisconnectedParts == FALSE) { 
    mygraph <- removeDisconnectedParts(mygraph)
    return(mygraph)
  }
  if (length(graph::nodes(mygraph)) == 0) {
    warning("Empty graph")
  }
  return(mygraph)
}
###############################################################################
#' This function is used internally by pathway2Graph to obtain physical entities participating in an interaction.
#'  
#' @author Nirupama Benis
#' @export
#' @import data.table
#' @examples

getParticipants <- function(pwComponentList, instance, biopaxlevel, splitComplexMolecules = FALSE, useIDasNodenames = TRUE) { 
  
  pwComponentList <- pwComponentList
  instance <- instance
  biopaxlevel <- biopaxlevel
  splitComplexMolecules <- splitComplexMolecules
  useIDasNodenames <- useIDasNodenames
  if (splitComplexMolecules & any(isOfClass(instance, "complex"))) {
    if (useIDasNodenames) {
      participants <- as.character(splitComplex(pwComponentList, instance$id[1], returnIDonly = TRUE, biopaxlevel = biopaxlevel))
    } else {
      participants <- as.character(splitComplex(pwComponentList, instance$id[1], biopaxlevel = biopaxlevel)$name)
    }
  } else {
    if (useIDasNodenames) {
      participants <- as.character(instance[1]$id)
    } else {
      participants <- getInstanceProperty(pwComponentList, instance[1]$id, biopaxlevel = biopaxlevel)
    }
  }
  return(participants)
}

###############################################################################
#' This function is used internally by pathway2Graph to remove the smaller disconnected parts of the pathway graph.
#'  
#' @author Nirupama Benis
#' @export
#' @import data.table
#' @examples

removeDisconnectedParts <- function(mygraph) { 
  mygraph <- mygraph
  if (!require(igraph)) {
    message(paste("Retrieving only the largest connected component needs the igraph library installed, albeit it cannot be found.", "\n"))
    return(NULL)
  }
  graphClusters <- clusters(igraph.from.graphNEL(mygraph))
  if(graphClusters$no > 1) { 
    for(j in 1:(graphClusters$no - 1)) { 
      smallestClust <- which.min(graphClusters$csize)
      smallestClustMembs <- names(graphClusters$membership)[graphClusters$membership == smallestClust]
      mygraph <- removeNode(smallestClustMembs, mygraph)
      graphClusters <- clusters(igraph.from.graphNEL(mygraph))
    }
  }
  return(mygraph)
}

###############################################################################
#' This function is used internally to remove the hash symbol.
#'  
#' @author Frank Kramer
#' @export
#' @import data.table
#' @examples

striphash <- function (x) {
  return(gsub("#", "", x))
}
