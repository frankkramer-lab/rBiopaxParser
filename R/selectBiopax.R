###############################################################################
#
# selectBiopax.R: 	This file contains the all functions related to selecting and retrieving information from a parsed Biopax model within R.
# author: Frank Kramer <dev@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################


#' Returns all list of all instances of a certain class, or all instances.
#' 
#' Returns a character vector of IDs of all instances of a certain class, or all instances if type=="all".
#' 
#' @param biopax A biopax model 
#' @param type string. Class of the instances
#' @returnType character vector
#' @return Returns a character vector containing all instances of the supplied class.
#' @author Frank Kramer
#' @export
getBiopaxInstancesList <- function(biopax, type="all") {
	if(tolower(type) == "all") {
		ret = unique(biopax$df[,c("instancetype", "instanceid")])
	} else {
		ret = unique(biopax$df[df$instancetype==stripns(type),c("instancetype", "instanceid")])
	}
	#ret[order(ret$instancetype,ret$instanceid),]
	ret
}

#' Returns all instances of a certain class.
#' 
#' Returns all instances of a certain class.
#' 
#' @param biopax A biopax model 
#' @param type string. Class of the instance
#' @returnType data.frame
#' @return Returns a data.frame containing all instances of the supplied class.
#' @author Frank Kramer
#' @export
getBiopaxInstancesByType <- function(biopax, type) {
	biopax$df[biopax$df$instancetype==stripns(type),]
}

#' Returns all instances with a certain ID.
#' 
#' Returns all instances with a certain ID.
#' 
#' @param biopax A biopax model 
#' @param id string. ID of the instance
#' @returnType data.frame
#' @return Returns a data.frame containing all instances with the supplied ID.
#' @author Frank Kramer
#' @export
getBiopaxInstancesByID <- function(biopax, id) {
	id = unique(striphash(id))
	biopax$df[biopax$df$instanceid %in% id,]
}

#' Returns all instances with a certain name.
#' 
#' Returns all instances with a certain name.
#' 
#' @param biopax A biopax model 
#' @param name string. Name of the instance
#' @returnType data.frame
#' @return Returns a data.frame containing all instances with the supplied name.
#' @author Frank Kramer
#' @export
getBiopaxInstancesByName <- function(biopax, name) {
	id = as.character(biopax$df[biopax$df$property=="NAME" & biopax$df$property_value == name,]$instanceid)
	id = unique(id)
	biopax$df[biopax$df$instanceid %in% id,]
}

#' Returns all ids of instances with a certain name.
#' 
#' Returns all ids of instances with a certain name.
#' 
#' @param biopax A biopax model 
#' @param name string. Name of the instance
#' @returnType character vector
#' @return Returns a vector containing all instance ids with the supplied name.
#' @author Frank Kramer
#' @export
getBiopaxIDsByName <- function(biopax, name) {
	id = as.character(biopax$df[biopax$df$property=="NAME" & biopax$df$property_value == name,]$instanceid)
	unique(id)
}

#' This function returns a list of all pathway ids.
#' 
#' This function returns a vector of all pathway ids.
#' 
#' @param biopax A biopax model 
#' @returnType character vector
#' @return Returns a character vector containing the names of all pathways.
#' @author Frank Kramer
#' @export
getPathwayList <- function(biopax) {
	#get pw ids, return a unique list
	biopax$df[tolower(biopax$df$instancetype) == "pathway" & tolower(biopax$df$property) == "name",c("instancetype", "instanceid", "property", "property_value")]
}

#' This function returns a data.frame containing all pathways.
#' 
#' This function returns a data.frame containing the all pathway instances.
#' 
#' @param biopax A biopax model 
#' @returnType data.frame
#' @return Returns a data.frame containing all pathway instances.
#' @author Frank Kramer
#' @export
getPathways <- function(biopax) {
	#list all pathways and return their data
	pathwayids = as.character(getPathwayList(biopax)$instanceid)
	getPathway(biopax, pathwayids)
}

#' This function returns a data.frame containing the pathway.
#' 
#' This function returns a data.frame containing the pathway instance pointed at by id.
#' 
#' @param biopax A biopax model
#' @param id string. ID of the pathway instance
#' @returnType data.frame
#' @return Returns a data.frame containing the pathway instance.
#' @author Frank Kramer
#' @export
getPathway <- function(biopax, id) {
	getBiopaxInstancesByID(biopax, c(getReferencedIDs(biopax,id),id) )
}

#' This function lists all pathway components of a given pathway.
#' 
#' This function returns a (unique) vector of pathway component IDs of the supplied pathway.
#' 
#' @param biopax A biopax model
#' @param id string. A pathway ID
#' @returnType character vector
#' @return Returns a character vector with IDs
#' @author Frank Kramer
#' @export
getPathwayComponentList <- function(biopax, id) {
	id = unique(striphash(id))
	#get pw component list
	pwcomp_list = as.character(biopax$df[tolower(biopax$df$property) == tolower("PATHWAY-COMPONENTS") & biopax$df$instanceid %in% id,"property_attr_value"])
	#strip # from front of id
	unique(striphash(pwcomp_list))
}

#' This function returns a data.frame containing all pathway components.
#' 
#' This function returns a data.frame containing all instances referenced in a pathways PATHWAY-COMPONENTS property.
#' 
#' @param biopax A biopax model 
#' @param id string. Pathway ID
#' @returnType data.frame
#' @return Returns a data.frame containing all pathway components.
#' @author Frank Kramer
#' @export
getPathwayComponents <- function(biopax, id) {
	id = unique(striphash(id))
	pwcomp_list = getPathwayComponentList(biopax,id)
	getBiopaxInstancesByID(biopax,pwcomp_list)
}

#' This function lists all components of a given complex.
#' 
#' This function returns a vector of all component IDs of the supplied complex.
#' 
#' @param biopax A biopax model
#' @param id string. A complex ID
#' @returnType character vector
#' @return Returns a character vector with IDs
#' @author Frank Kramer
#' @export
getComplexComponentList <- function(biopax, id) {
	unlist(getReferencedIDs(biopax, id, recursive=FALSE, onlyFollowProperties=c("COMPONENTS")))
}

#' This functions splits up a complex into its components.
#' 
#' This function looks up the supplied Complex ID and returns the names of all its components. 
#' 
#' @param biopax A biopax model 
#' @param complexid string ID of an complex
#' @param recursive logical
#' @returnType character vector
#' @return Returns a character vector with the names of all subcomponents.
#' @author Frank Kramer
#' @export
splitComplex <- function(biopax, complexid, recursive=TRUE) {
	# complexes can contain entries via "bp:COMPONENTS" -> physicalentityparticipant "bp:PHYSICAL-ENTITY" -> physicalentity
	referenced = getReferencedInstances(biopax, complexid, recursive=recursive, onlyFollowProperties=c("COMPONENTS","PHYSICAL-ENTITY"))
	
	sel = tolower(referenced$instancetype) %in% tolower(c("dna","rna","protein","smallMolecule"))
	unique(referenced[referenced$property=="NAME" & sel,c("instanceid","property_value")])
	
#	#every component is a physicalentityparticipant
#	peps = getBioPaxInstancesByID(df, instance[instance$property=="bp:COMPONENTS","property_attr_value"])
#	#each physicalentityparticipant has exactly 1 physicalentity, get that and get the name!
#	pes = getBioPaxInstancesByID(df, peps[peps$property=="bp:PHYSICAL-ENTITY","property_attr_value"])
#	
#	if(recursive) {
#		sel_complex = isOfClass(pes,class="bp:complex")
#		ret = unique(pes[pes$property=="bp:NAME" & !sel_complex,c("instanceid","property_value")])
#		while(any(sel_complex)) {
#			for( subid in unique(pes$instanceid[sel_complex]) ) {
#				
#			}
#		}
#	} else {
#		ret = unique(pes[pes$property=="bp:NAME",c("instanceid","property_value")])		
#	}
#	ret
}


#' This function returns a vector of ids of all referenced instances of the supplied instance.
#' 
#' This function takes an id and a biopax model as input. The id of every instance that is referenced is returned.
#' If recursive == TRUE this function recurses through all referenced IDs of the referenced instances and so on.
#' "onlyFollowProperties" limits the recursivness to only certain properties, for example follow only complexes or physicalEntities.
#' 
#' @param biopax A biopax model 
#' @param id string. ID of the instance
#' @param recursive logical
#' @param onlyFollowProperties character vector
#' @returnType character vector
#' @return Returns a character vector of IDs referenced by the supplied id in the supplied biopax model.
#' @author Frank Kramer
#' @export
getReferencedIDs <- function(biopax, id, recursive=TRUE, onlyFollowProperties=c()) {
	id = unique(striphash(id))
	referencedIDs = list()
	
	#every ref in instances of id
	if(length(onlyFollowProperties) > 0) {
		newIDs = biopax$df[biopax$df$instanceid %in% id & biopax$df$property_attr == "rdf:resource" & tolower(biopax$df$property) %in% tolower(onlyFollowProperties),"property_attr_value"]	
	} else {
		newIDs = biopax$df[biopax$df$instanceid %in% id & biopax$df$property_attr == "rdf:resource" ,"property_attr_value"]	
	}
	newIDs = unique(striphash(newIDs))
	newIDs = newIDs[!(newIDs %in% id)]
	referencedIDs = c(referencedIDs,newIDs)
	
	if(recursive) {
		while(length(newIDs)>0) {
			if(length(onlyFollowProperties) > 0) {
				newIDs = biopax$df[biopax$df$instanceid %in% newIDs & biopax$df$property_attr == "rdf:resource" & tolower(biopax$df$property) %in% tolower(onlyFollowProperties),"property_attr_value"]	
			} else {
				newIDs = biopax$df[biopax$df$instanceid %in% newIDs & biopax$df$property_attr == "rdf:resource" ,"property_attr_value"]	
			}
			newIDs = unique(striphash(newIDs))
			newIDs = newIDs[!(newIDs %in% c(referencedIDs,id))]
			referencedIDs = c(referencedIDs,newIDs)
		}
	}
	referencedIDs
}

#' This function returns a data.frame containing all referenced instances of the supplied instance.
#' 
#' This function takes an id and a biopax model as input. A data.frame containing all instances referenced by the supplied instance is returned.
#' If recursive == TRUE this function recurses through all referenced IDs of the referenced instances and so on. 
#' "onlyFollowProperties" limits the recursivness to only certain properties, for example follow only complexes or physicalEntities.
#' 
#' @param biopax A biopax model
#' @param id string
#' @param recursive logical
#' @param onlyFollowProperties character vector
#' @returnType data.frame
#' @return Returns a data.frame containing all referenced instances.
#' @author Frank Kramer
#' @export
getReferencedInstances <- function(biopax, id, recursive=TRUE, onlyFollowProperties=c()) {
	referencedIDs = getReferencedIDs(biopax, id, recursive, onlyFollowProperties)
	getBiopaxInstancesByID(biopax, referencedIDs)
}

#' This function returns a vector of ids of all instances that reference the supplied instanceid.
#' 
#' This function takes an id and a biopax model as input. The id of every instance that references the supplied id is returned.
#' If recursive == TRUE this function recurses through all referencing IDs of the referencing instances and so on.
#' "onlyFollowProperties" limits the recursivness to only certain properties, for example follow only complexes or physicalEntities.
#' 
#' @param biopax A biopax model 
#' @param id string. ID of the instance
#' @param recursive logical
#' @param onlyFollowProperties character vector
#' @returnType character vector
#' @return Returns a character vector of IDs referencing the supplied id in the supplied biopax model.
#' @author Frank Kramer
#' @export
getReferencingIDs <- function(biopax, id, recursive=TRUE, onlyFollowProperties=c()) {
	id = unique(striphash(id))
	id = addhash(id)
	referencingIDs = list()
	#####################################continue here! how to traverse backwards through biopax networks??
	#every ref in instances of id
	if(length(onlyFollowProperties) > 0) {
		newIDs = biopax$df[biopax$df$property_attr_value %in% id & biopax$df$property_attr == "rdf:resource" & tolower(biopax$df$property) %in% tolower(onlyFollowProperties),"instanceid"]	
	} else {
		newIDs = biopax$df[biopax$df$property_attr_value %in% id & biopax$df$property_attr == "rdf:resource" ,"instanceid"]	
	}
	newIDs = unique(addhash(newIDs))
	newIDs = newIDs[!(newIDs %in% id)]
	referencingIDs = c(referencingIDs,newIDs)
	
	if(recursive) {
		while(length(newIDs)>0) {
			if(length(onlyFollowProperties) > 0) {
				newIDs = biopax$df[biopax$df$property_attr_value %in% newIDs & biopax$df$property_attr == "rdf:resource" & tolower(biopax$df$property) %in% tolower(onlyFollowProperties),"instanceid"]	
			} else {
				newIDs = biopax$df[biopax$df$property_attr_value %in% newIDs & biopax$df$property_attr == "rdf:resource" ,"instanceid"]	
			}
			newIDs = unique(striphash(newIDs))
			newIDs = newIDs[!(newIDs %in% c(referencingIDs,id))]
			referencingIDs = c(referencingIDs,newIDs)
		}
	}
	striphash(referencingIDs)
}

#' This function returns a data.frame containing all instances referencing the supplied instance.
#' 
#' This function takes an id and a biopax model as input. A data.frame containing all instances referencing the supplied instance is returned.
#' If recursive == TRUE this function recurses through all referencing IDs of the referencing instances and so on. 
#' "onlyFollowProperties" limits the recursivness to only certain properties, for example follow only complexes or physicalEntities.
#' 
#' @param biopax A biopax model
#' @param id string
#' @param recursive logical
#' @param onlyFollowProperties character vector
#' @returnType data.frame
#' @return Returns a data.frame containing all referencing instances.
#' @author Frank Kramer
#' @export
getReferencingInstances <- function(biopax, id, recursive=TRUE, onlyFollowProperties=c()) {
	referencedIDs = getReferencedIDs(biopax, id, recursive, onlyFollowProperties)
	getBiopaxInstancesByID(biopax, referencedIDs)
}

#' This function returns the class name of the instance.
#' 
#' This function returns the class name of the instance.
#' 
#' @param biopax A biopax model
#' @param instanceid string
#' @returnType string
#' @return Returns the class name of the biopax instance.
#' @author fkramer
#' @export
getInstanceClass <- function(biopax, instanceid) {
	instanceid = striphash(instanceid)
	as.character(biopax$df[biopax$df$instanceid == instanceid,"instancetype"][1])
}

#' This function returns all properties of the specified type for an instance.
#' 
#' This function returns all properties of the specified type for an instance.
#' 
#' @param biopax A biopax model
#' @param instanceid string
#' @param property string. Attention: All properties in Biopax Level 2 are all upper case. 
#' @returnType character vector
#' @return Returns a character vector with all properties of the selected type for this instance. Returns NULL if no property of this type is found.
#' @author fkramer
#' @export
getInstanceProperty <- function(biopax, instanceid, property="NAME") {
	instanceid = striphash(instanceid)
	#get class of the instance
	instancetype = getInstanceClass(biopax,instanceid)
	#get the properties this instance has
	classproperties = getClassProperties(instancetype)
	#if its a reference return the property_attr_value (its the referenced ids) or return property_value otherwise
	property_type = unlist(classproperties[classproperties$property == property,]$property_type)[1]
	
	#value
	if(length(property_type)>0) {
		if(property_type %in% c("http://www.w3.org/2001/XMLSchema#string","http://www.w3.org/2001/XMLSchema#double","http://www.w3.org/2001/XMLSchema#float", "http://www.w3.org/2001/XMLSchema#integer" )) {
			as.character(biopax$df[biopax$df$instanceid == instanceid & biopax$df$property == property,"property_value"])
		} else { #reference
			as.character(biopax$df[biopax$df$instanceid == instanceid & biopax$df$property == property,"property_attr_value"])
		}
	} else {
		return(NULL)
	}
}

#' This function resolves physicalEntityParticipantIDs to their corresponding physicalEntityIDs
#' 
#' This function resolves physicalEntityParticipantIDs to their corresponding physicalEntityIDs. Every physicalEntityParticipant corresponds exactly to one physicalEntity.
#' 
#' @param biopax A biopax model
#' @param physicalEntityId string. IDs of physicalEntityParticipants to be resolved
#' @returnType character vector 
#' @return Returns ids of physicalEntity corresponding to the specified  physicalEntityParticipantIDs
#' @author fkramer
internal_resolvePhysicalEntityParticipant <- function(biopax, physicalEntityId) {
	unlist(getReferencedIDs(biopax, physicalEntityId, recursive=FALSE, onlyFollowProperties=c("PHYSICAL-ENTITY")))
}

#' This function returns the neighborhood of a physicalEntity
#' 
#' This function searches the supplied biopax for interactions that are connected to the molecule or within 'depth' number of steps from it.
#' 
#' @param biopax A biopax model
#' @param instanceid string. ID of a physicalEntity (dna, rna, protein, complex, smallMolecule)
#' @param depth integer. Search depth, this specifies how far out from the specified molecule the neighborhood should be streched.
#' @param onlyInPathways character vector of pathway IDs. Search only in these pathways for neighbors. 
#' @returnType character vector 
#' @return Returns ids of interactions within 'depth' number of steps of the specified physicalEntity 
#' @author fkramer
#' @export
getNeighborhood <- function(biopax, instanceid, depth=1, onlyInPathways=c()) {
	interactionlist = c()
	moleculelist = instanceid
	for(i in 1:depth) {
		# add all interaction that includes at least one of the molecules to the interactionlist
		#biopax$df[biopax$df$instanceid == instanceid & biopax$df$property == property,"property_attr_value"]
		getReferencedInstances(biopax, instanceid, recursive=TRUE, onlyFollowProperties=c("COMPONENTS","PHYSICAL-ENTITY"))
	}
}

#' This function returns the annotations of the supplied IDs
#' 
#' This function returns the annotations of the supplied IDs in a data.frame.
#' 
#' @param biopax A biopax model
#' @param id vector of strings. IDs of instances to get annotations
#' @param splitComplexes logical. If TRUE complexes are split up into their components and the annotation of the components is added.
#' @param followPhysicalEntityParticipants logical. If TRUE physicalEntityParticipants are resolved to their corresponding physicalEntities and their annotation is added. 
#' @returnType data.frame
#' @return Returns data.frame with annotations 
#' @author fkramer
#' @export
getXrefAnnotations <- function(biopax, id, splitComplexes=FALSE, followPhysicalEntityParticipants=TRUE) {
	id = c(striphash(id))
	#annotations = data.frame(type=NA, id=NA, name=NA, annotation_type=NA, annotation_id=NA, annotation=NA, stringsAsFactors=FALSE)
	annotations = matrix(nrow=0, ncol=6)
	colnames(annotations) = c("type", "id","name","annotation_type", "annotation_id", "annotation")
	
	for(i in 1:length(id)) {
		# if its a complex AND we're supposed to split it:
		if(getInstanceClass(biopax,id[i]) == "complex" & splitComplexes) {
			referenced = getReferencedInstances(biopax, id[i], recursive=TRUE, onlyFollowProperties=c("COMPONENTS","PHYSICAL-ENTITY"))
			sel = tolower(referenced$instancetype) %in% tolower(c("dna","rna","protein","smallMolecule"))
			referenced = as.character(unique(referenced[sel & !(referenced$instanceid %in% id),"instanceid"]))
			annotations = rbind(annotations, getXrefAnnotations(biopax,referenced, splitComplexes=FALSE, followPhysicalEntityParticipants=FALSE))
		} else if(getInstanceClass(biopax,id[i]) == "physicalEntityParticipant" & followPhysicalEntityParticipants) {
			# if its a physicalEntityParticipant
			peID = internal_resolvePhysicalEntityParticipant(biopax, id[i])
			if(!(peID %in% id)) {
				annotations = rbind(annotations, getXrefAnnotations(biopax, peID, splitComplexes=splitComplexes, followPhysicalEntityParticipants=TRUE))
			}
		} else {
			# for any other class do this
			type = getInstanceClass(biopax, id[i])
			name = getInstanceProperty(biopax, id[i], property="NAME")[1]
			xrefs = getInstanceProperty(biopax, id[i], property="XREF")
			for(xref in xrefs) {
				annotations = rbind(annotations, c(type, id[i], name, getInstanceClass(biopax, xref), striphash(xref),
								paste(getInstanceProperty(biopax, xref, property="DB"),	":",
										getInstanceProperty(biopax, xref, property="ID"), sep="")
								))
			}
		}
	}
	annotations
}


########## NO SUPPORT FOR MANIPULATING XML NODES DIRECTLY AT THE MOMENT!!!
#
#getBioPaxInstances.owl <- function(owl, type) {
#	if( !(tolower(substring(type,1,3)) == "bp:") ) {
#		type = paste("bp:",type,sep="")
#	}
#	getNodeSet(owl$biopaxxml,paste("/rdf:RDF/",type,sep=""))
#}
#
#getBioPaxInstance.owl <- function(owl, id) {
#	getNodeSet(owl$biopaxxml,paste("/rdf:RDF//*[@*='",id,"']",sep=""))		
#}
#
#getPathway.owl <- function(owl, id) {
#	
#}
#

