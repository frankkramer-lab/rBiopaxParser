###############################################################################
#
# modifyBiopax.R: 	This file contains the all functions related to modifying a parsed Biopax model within R.
# author: Frank Kramer <mail@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################


#' This function adds new instances to an existing biopax model.
#' 
#' This function adds new instances (supplied as a compatible data.frame) to an existing biopax model.
#' 
#' @param biopax A biopax model 
#' @param newInstancesDF data.frame. Compatible with internal Biopax Level 2 implementation.
#' @returnType A biopax model
#' @return Returns the supplied biopax model with the new instances added.
#' @author Frank Kramer
#' @export
addBiopaxInstances <- function(biopax,newInstancesDF) {
	biopax$df = rbind(biopax$df, newInstancesDF)
	biopax
}

#' This function adds a new instance an existing biopax model.
#' 
#' This function adds a new instance an existing biopax model.
#' "properties" is a named list of vectors, with the vector name as the name of the property and every entry of the vector a property value.
#' Please note: case sensitivity! In Biopax Level 2 all properties are written in all capital letters.
#' 
#' @param biopax A biopax model 
#' @param instancetype string. Class name
#' @param instanceid string. ID of the instance
#' @param properties named list of properties.
#' @param verbose logical. Be verbose about what was added.
#' @returnType A biopax model
#' @return Returns the supplied biopax model with the new instance added.
#' @author Frank Kramer
#' @export
addBiopaxInstance <- function(biopax, instancetype, instanceid, properties=list(NAME=c()), verbose=TRUE) {
	propertyDF = internal_propertyListToDF(instancetype, instanceid, properties, namespace_rdf=biopax$ns_rdf)
	biopax$df = rbind(biopax$df,propertyDF)
	if(verbose) cat(paste("Added", instancetype, "with ID:", instanceid,"\n"))
	biopax
}

#' This function adds new properties to an existing biopax instance.
#' 
#' This function adds new properties to an existing biopax instance.
#' 
#' @param biopax A biopax model
#' @param instanceid string. ID of the instance
#' @param properties named list of properties.
#' @returnType A biopax model
#' @return Returns the supplied biopax model with new properties added to this instance.
#' @author Frank Kramer
#' @export
addPropertiesToBiopaxInstance <- function(biopax,instanceid, properties) {
	instancetype = getInstanceClass(biopax,instanceid)
	propertyDF = internal_propertyListToDF(instancetype, instanceid, properties, namespace_rdf=biopax$ns_rdf)
	biopax$df = rbind(biopax$df,propertyDF)
	biopax
}

#' Internal function to build a data.frame from the list of properties for a new instance
#' 
#' Internal function to build a data.frame from the list of properties for a new instance
#' 
#' @param instancetype string. Class name
#' @param instanceid string. ID of the instance
#' @param properties named list of properties.
#' @param namespace_rdf string. This defines the rdf namespace to use.
#' @returnType data.frame
#' @return Returns a data.frame with the new properties for the given instance
#' @author Frank Kramer
#  #' @export
internal_propertyListToDF <- function(instancetype, instanceid, properties, namespace_rdf="rdf") {
	
	ret = matrix(data=NA,nrow = 1000, ncol = 6)
	ret_colnames = c("instancetype","instanceid","property","property_attr","property_attr_value","property_value")
	colnames(ret) = ret_colnames
	
	#every property of an instance is represented as 1 entry in the df
	rowcount = 0
	names(properties) = toupper(names(properties))
	
	#
	classproperties = getClassProperties(instancetype)
	
	for(n in names(properties)) {
		
		# find out if this is a value or a reference
		property_type = unlist(classproperties[classproperties$property == n,]$property_type)[1]
		
		for(p in properties[[n]]) {
			rowcount = rowcount + 1
			
			#value
			if(property_type %in% c("http://www.w3.org/2001/XMLSchema#string","http://www.w3.org/2001/XMLSchema#double","http://www.w3.org/2001/XMLSchema#float", "http://www.w3.org/2001/XMLSchema#integer" )) {
				ret[rowcount,] = c(
						instancetype, instanceid,
						n, #property				
						paste(namespace_rdf,":datatype",sep=""), #property_attr
						property_type, #property_attr_value
						p #property_value
				)
			} else { #reference
				ret[rowcount,] = c(
						instancetype, instanceid,
						n, #property				
						paste(namespace_rdf,":resource",sep=""), #property_attr
						paste("#",p,sep=""), #property_attr_value
						"" #property_value
				)
			}
		}
	}
	matrix( ret[!is.na(ret[,1]),], ncol=6, dimnames=list(list(),ret_colnames))
	
}


#' This function adds a new pathway to the biopax model.
#' 
#' This function adds a new pathway + its PATHWAY-COMPONENTS (references to interaction/pathways/pathwaySteps)
#' 
#' @param biopax A biopax model 
#' @param NAME string. Name of the pathway
#' @param PATHWAY_COMPONENTS character vector. IDs of the pathway components. This must be IDs of instances of type interaction/pathway/pathwayStep (or their subclasses).
#' @param ID string. ID for the pathway. If NULL a new ID is generated with prefix "pathway".
#' @param ORGANISM string. Organism property of the pathway. optional.
#' @param COMMENT string. An optional comment 
#' @returnType biopax
#' @return Returns the biopax model with the added pathway.
#' @author fkramer
#' @export
addPathway <- function(biopax, NAME, PATHWAY_COMPONENTS=c(), ID=NULL, ORGANISM=NULL, COMMENT=NULL) {
	
	if(is.null(ID)) ID = generateNewUniqueID(biopax, id="pathway")
	if( ID %in% biopax$df$instanceid ) ID = generateNewUniqueID(biopax, id=ID)
	
	properties = list(NAME=c(NAME),'PATHWAY-COMPONENTS'=PATHWAY_COMPONENTS)
	if(!is.null(ORGANISM)) properties['ORGANISM']=c(ORGANISM)
	if(!is.null(COMMENT)) properties['COMMENT']=c(COMMENT)
	
	biopax = addBiopaxInstance(biopax, instancetype="pathway", instanceid=ID, properties=properties)
	biopax
}

#' This function adds pathway components to an existing pathway
#' 
#' This function adds pathway components to an existing pathway.
#' Property PATHWAY-COMPONENTS  are references to IDs of interaction/pathways/pathwaySteps (or subclasses of those)
#' 
#' @param biopax A biopax model
#' @param ID string. ID for the pathway
#' @param PATHWAY_COMPONENTS character vector. IDs of the pathway components. This must be IDs of instances of type interaction/pathway/pathwayStep (or their subclasses).
#' @returnType biopax
#' @return Returns the biopax model with the pathway components added to the pathway
#' @author fkramer
#' @export
addPathwayComponents <- function(biopax, ID, PATHWAY_COMPONENTS=c()) {
	properties = list('PATHWAY-COMPONENTS'=PATHWAY_COMPONENTS)
	biopax = addPropertiesToBiopaxInstance(biopax, instanceid=ID, properties=properties)
	biopax
}

#' This function adds a new control to the biopax model.
#' 
#' This function adds a new interaction of class control to the biopax model. This is a convenience function to add controls,
#' internally the function addBiopaxInstance is called with properties CONTROL-TYPE, CONTROLLER and CONTROLLED set.
#' 
#' @param biopax A biopax model
#' @param CONTROL_TYPE string. Specifies wether this is an activating or inhibiting control.
#' @param CONTROLLER string. ID of the physicalEntityParticipant instance that is the controller of this interaction.
#' @param CONTROLLED vector of strings. IDs of the interaction and/or pathway instances that are being controlled.
#' @param ID string. ID for the control. If NULL a new ID is generated with prefix "control". 
#' @returnType biopax
#' @return Returns the biopax model with the added pathway.
#' @author fkramer
#' @export
addControl <- function(biopax, CONTROL_TYPE=c("ACTIVATION","INHIBITION"), CONTROLLER="", CONTROLLED=c(), ID=NULL) {
	
	if(is.null(ID)) ID = generateNewUniqueID(biopax, id="control")
	if( ID %in% biopax$df$instanceid ) ID = generateNewUniqueID(biopax, id=ID)
	
	properties = list('CONTROL-TYPE'=c(CONTROL_TYPE[1]), CONTROLLER=c(CONTROLLER), CONTROLLED=CONTROLLED)
	
	biopax = addBiopaxInstance(biopax, instancetype="control", instanceid=ID, properties=properties)
	biopax
}

#' This function adds a new biochemical reaction to the biopax model.
#' 
#' This function adds a new biochemical reaction of class biochemicalReaction to the biopax model. This is a convenience function,
#' internally the function addBiopaxInstance is called with properties LEFT and RIGHT set.
#' 
#' @param biopax A biopax model
#' @param LEFT vector of strings. IDs of the physicalEntityParticipant instances that are on the left side of this reaction.
#' @param RIGHT vector of strings. IDs of the physicalEntityParticipant instances that are on the right side of this reaction.
#' @param ID string. ID for the control. If NULL a new ID is generated with prefix "biochemicalReaction". 
#' @returnType biopax
#' @return Returns the biopax model with the added pathway.
#' @author fkramer
#' @export
addBiochemicalReaction <- function(biopax, LEFT=c(), RIGHT=c(), ID=NULL) {
	
	if(is.null(ID)) ID = generateNewUniqueID(biopax, id="biochemicalReaction")
	if( ID %in% biopax$df$instanceid ) ID = generateNewUniqueID(biopax, id=ID)
	
	properties = list(LEFT=c(LEFT), RIGHT=c(RIGHT))
	
	biopax = addBiopaxInstance(biopax, instancetype="biochemicalReaction", instanceid=ID, properties=properties)
	biopax
}

#' This function adds a new physical entity.
#' 
#' This function adds a new physical entity of chosen class to the biopax model. This is a convenience function to add physical entities,
#' internally the function addBiopaxInstance is called with properties NAME and ORGANISM set.
#' 
#' @param biopax A biopax model
#' @param class string. Class of the physical entity to add, choose from c("dna","rna","protein","smallMolecule","complex").  
#' @param NAME string. Name of the new physical entity
#' @param ID string. ID for the physical entity. If NULL a new ID is generated with prefix "physicalEntity".
#' @param ORGANISM string. Organism property of the molecule. optional.
#' @param COMMENT string. An optional comment 
#' @returnType biopax
#' @return Returns the biopax model with the added physical entity.
#' @author fkramer
#' @export
addPhysicalEntity <- function(biopax, class=c("dna","rna","protein","smallMolecule","complex")[1], NAME, ID=NULL, ORGANISM=NULL, COMMENT=NULL) {
	
	if(is.null(ID)) ID = generateNewUniqueID(biopax, id="physicalEntity")
	if( ID %in% biopax$df$instanceid ) ID = generateNewUniqueID(biopax, id=ID)
	
	properties = list(NAME=c(NAME))
	if(!is.null(ORGANISM)) properties['ORGANISM']=c(ORGANISM)
	if(!is.null(COMMENT)) properties['COMMENT']=c(COMMENT)
	
	biopax = addBiopaxInstance(biopax, instancetype=class, instanceid=ID, properties=properties)
	biopax
}

#' This function adds a new physical entity participant.
#' 
#' This function adds a new physical entity participant instance, which is a placeholder for physicalEntity class instances in interactions. 
#' This is a convenience function to add physicalEntityParticipant instances, internally the function addBiopaxInstance is called.
#' 
#' @param biopax A biopax model
#' @param referencedPhysicalEntityID string. ID the new physicalEntity instance to reference here.
#' @param ID string. ID for the physical entity participant. If NULL a new ID is generated with prefix "physicalEntityParticipant".
#' @returnType biopax
#' @return Returns the biopax model with the added physicalEntityParticipant.
#' @author fkramer
#' @export
addPhysicalEntityParticipant <- function(biopax, referencedPhysicalEntityID, ID=NULL) {
	
	if(is.null(ID)) ID = generateNewUniqueID(biopax, id="physicalEntityParticipant")
	if( ID %in% biopax$df$instanceid ) ID = generateNewUniqueID(biopax, id=ID)
	
	properties = list('PHYSICAL-ENTITY'=c(referencedPhysicalEntityID))
	
	biopax = addBiopaxInstance(biopax, instancetype="physicalEntityParticipant", instanceid=ID, properties=properties)
	biopax
}

#' This function generates a new unique id for a biopax model
#' 
#' This function generates a new unique id for a biopax model. Pass it an startin g point like "pathway" or "protein" to get a niceer looking id. 
#' 
#' @param biopax A biopax model
#' @param id string. This is used as a prefix for the id.
#' @returnType string
#' @return Returns an unused unique ID.
#' @author fkramer
#' @export
generateNewUniqueID <- function(biopax, id="") {
	if(id=="") id="id"
		
	for(i in 1:500000) {
		myid = paste(id,i,sep="_")
		if( !(myid %in% biopax$df$instanceid) ) {
			return(myid)
		}
	}
	stop("Error: generateNewUniqueID ran out of new ideas for ids :-(")
}


#' This function merges two given pathways 
#' 
#' This function merges two given pathways and appends it to the supplied biopax model. The user has to specify a new name for the pathways and can
#' supply ID, ORGANISM and COMMENT properties for the new pathway. If no ID is supplied, a new unique ID is generated. If no organism property is supplied
#' the organism property of the first pathway is re-used. If ORGANISM is NULL the property is not set. Optionally a comment can be added to the pathway.   
#' 
#' @param biopax A biopax model 
#' @param pwid1 string. ID of first pathway to merge
#' @param pwid2 string. ID of second pathway to merge
#' @param NAME string. Name of the new merged pathway
#' @param ID string. ID for the pathway. If NULL a new ID is generated with prefix "pathway".
#' @param ORGANISM string. Organism property of the pathway. By default uses the same organism as the first supplied pathway. If NULL no organism property is set.
#' @param COMMENT string. An optional comment 
#' @returnType A biopax model
#' @return A biopax model with the merged pathway added.
#' @author fkramer
#' @export
mergePathways <- function(biopax, pwid1, pwid2, NAME, ID=NULL, ORGANISM="", COMMENT=NULL) {
	
	if(is.null(ID)) ID = generateNewUniqueID(biopax, id="pathway")
	if( ID %in% biopax$df$instanceid ) ID = generateNewUniqueID(biopax, id=ID)
	
	PATHWAY_COMPONENTS = unique(c(getPathwayComponentList(biopax,pwid1),getPathwayComponentList(biopax,pwid2)))
	
	properties = list(NAME=c(NAME),'PATHWAY-COMPONENTS'=PATHWAY_COMPONENTS)
	if(!is.null(ORGANISM)) if(ORGANISM=="") ORGANISM=getInstanceProperty(biopax, instanceid=pwid1, property="ORGANISM")
	if(!is.null(ORGANISM)) if(length(ORGANISM)>0) properties['ORGANISM']=c(ORGANISM)
	if(!is.null(COMMENT)) properties['COMMENT']=c(COMMENT)
	
	biopax = addBiopaxInstance(biopax, instancetype="pathway", instanceid=ID, properties=properties)
	biopax
	
}

#' This function removes an instance
#' 
#' This function removes an instance from an existing biopax model.
#' 
#' @param biopax A biopax model
#' @param instanceid string. ID of the instance
#' @returnType A biopax model
#' @return Returns the supplied biopax model with the instance removed from it.
#' @author Frank Kramer
#' @export
removeInstance <- function(biopax, instanceid) {
	biopax$df = biopax$df[biopax$df$instanceid != instanceid, ]
	biopax
}

#' This function removes a property
#' 
#' This function removes a property fram an existing biopax instance.
#' 
#' @param biopax A biopax model
#' @param instanceid string. ID of the instance
#' @param properties character vector. listing the properties to remove.
#' @returnType A biopax model
#' @return Returns the supplied biopax model with properties removed from this instance.
#' @author Frank Kramer
#' @export
removeProperties <- function(biopax, instanceid, properties) {
	biopax$df = biopax$df[biopax$df$instanceid != instanceid | !(biopax$df$property %in% properties), ]
	biopax
}
