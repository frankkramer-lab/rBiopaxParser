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


#' Returns all instances that conform to the selection criteria.
#' 
#' Returns all instances that conform to the selection criteria. This function returns a subset of the internal data.frame of the biopax object.
#' Selection criteria are wether instances belong to a certain class or have the specified id, property or name. Setting a criteria to NULL ignores this criteria.
#' If returnValues is set to FALSE only the selector (a logical vector with length of the internal data.frame) is returned, otherwise the selected data is returned.
#' If includeSubClasses is set to TRUE the class criteria is broadened to include all classes that inherit from the given class, e.g. if class="control" and includeSubClasses=TRUE the function will select catalyses and modulations too, since they are a subclass of class control. 
#' If includeReferencedInstances is set to TRUE all instances that are being referenced by the selected instances are being selected too. The parameter works recursively, this means for example that a selected pathway and all it's interactions, complexes, molecules and annotations are returned if this parameter is set to true. This parameter is especially helpful if you want to migrate or merge knowledge from different data bases.
#' 
#' @param biopax A biopax model
#' @param class string. Class of the instances to select
#' @param id string. ID of the instances to select
#' @param property string. Return only this property of the instances
#' @param name string. Name of the instances to select
#' @param returnValues logical. If returnValues is set to FALSE only the selector (a logical vector with length of the internal data.frame) is returned, otherwise the selected data is returned
#' @param includeSubClasses logical. If includeSubClasses is set to TRUE the class criteria is broadened to include all classes that inherit from the given class
#' @param includeReferencedInstances logical. If includeReferencedInstances is set to TRUE all instances that are being referenced by the selected instances are being selected too
#' @return Returns a data.frame containing all instances conforming to the given selection criteria if returnValues=TRUE, only the selector for the internal data.frame otherwise.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  # select the subset of the internal data.frame that belongs to class "protein"
#'  selectInstances(biopax, class="protein")
#'  # select the subset of the internal data.frame that belongs to class "interaction"
#'  selectInstances(biopax, class="interaction")
#'  # select the subset of the internal data.frame that belongs to class "interaction" or any of its sub classes, like control, catalysis etc.
#'  selectInstances(biopax, class="interaction", includeSubClasses=TRUE)
#'  # select the subset of the internal data.frame that belongs to class "pathway" AND is a "NAME" property
#'  selectInstances(biopax, class="pathway", property="NAME")
selectInstances <- function (biopax, class=NULL, id=NULL, property=NULL, name=NULL, returnValues=TRUE, includeSubClasses=FALSE, includeReferencedInstances=FALSE) {
	sel = rep.int(TRUE,dim(biopax$df)[1])
	
	if(!is.null(id)) {
		id = unique(striphash(id))
		sel = sel & (biopax$df$id %in% id)
	}
	
	if(!is.null(class)) {
		if(includeSubClasses) {
			class = unique(c(class,getSubClasses(class, biopax$biopaxlevel)))
		}
		sel = sel & (tolower(biopax$df$class) %in% tolower(stripns(class)))
	}
	
	if(!is.null(property)) {
		sel = sel & (tolower(biopax$df$property) %in% tolower(property))
	}
	
	if(!is.null(name)) {
		ids = as.character(biopax$df[tolower(biopax$df$property)=="name" & biopax$df$property_value %in% name,]$id)
		ids = unique(ids)
		sel = sel & (biopax$df$id %in% ids)
	}
	
	if(includeReferencedInstances) {
		ids = as.character(biopax$df[sel,]$id)
		ids = unique(ids)
		ids = getReferencedIDs(biopax, ids)
		sel = sel & (biopax$df$id %in% ids)
	}	
	
	if(returnValues) {
		biopax$df[sel,]
	} else {
		sel
	}
	
}


#' Lists all instances that conform to the selection criteria.
#' 
#' Lists all instances that conform to the selection criteria. In contrast to selectInstances this function returns an easier to read list.
#' This function returns an ordered data.frame of class, id and name of the instances.
#' Selection criteria are wether instances belong to a certain class or have the specified id or name. Setting a criteria to NULL ignores this criteria.
#' If includeSubClasses is set to TRUE the class criteria is broadened to include all classes that inherit from the given class, e.g. if class="control" and includeSubClasses=TRUE the function will select catalyses and modulations too, since they are a subclass of class control. 
#' 
#' @param biopax A biopax model
#' @param class string. Class of the instances to select
#' @param id string. ID of the instances to select
#' @param name string. Name of the instances to select
#' @param includeSubClasses logical. If includeSubClasses is set to TRUE the class criteria is broadened to include all classes that inherit from the given class
#' @return Returns a data.frame containing all instances conforming to the given selection criteria if returnValues=TRUE, only the selector for the internal data.frame otherwise.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  # list all instances of class "protein"
#'  listInstances(biopax, class="protein")
#'  # list all instances of class "pathway"
#'  listInstances(biopax, class="pathway")
#'  # list all interaction including all subclasses of interactions
#'  listInstances(biopax, class="interaction", includeSubClasses=TRUE)
listInstances <- function (biopax, class=NULL, id=NULL, name=NULL, includeSubClasses=FALSE) {
	sel = rep.int(TRUE,dim(biopax$df)[1])

	if(!is.null(id)) {
		id = unique(striphash(id))
		sel = sel & (biopax$df$id %in% id)
	}
	
	if(!is.null(class)) {
		if(includeSubClasses) {
			class = unique(c(class,getSubClasses(class, biopax$biopaxlevel)))
		}
		sel = sel & (tolower(biopax$df$class) %in% tolower(stripns(class)))
	}
	
	if(!is.null(name)) {
		ids = as.character(biopax$df[tolower(biopax$df$property)=="name" & biopax$df$property_value %in% name,]$id)
		ids = unique(ids)
		sel = sel & (biopax$df$id %in% ids)
	}
	
	ret = unfactorize(unique(biopax$df[sel,c("class","id")]))
	if(dim(ret)[1]==0) return(NULL)
	ret$name = vector(mode="character", length=length(ret$id))
	colnames(ret) = c("class","id","name")
	#get name property for these instances. this slows things down but names seem so important ;-)
	subbiopax = biopax
	subbiopax$df = biopax$df[sel,]
	for(i in 1:length(ret$id)) {
		n = getInstanceProperty(subbiopax, ret$id[i])
		if(!is.null(n) && !length(n)==0) ret$name[i] = n
	}	
	
	ret
}


#' This function returns a list of all pathway ids.
#' 
#' This function returns a vector of all pathway ids.
#' 
#' @param biopax A biopax model 
#' @return Returns a character vector containing the names of all pathways.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  listPathways(biopax)
listPathways <- function(biopax) {
	listInstances(biopax, class="pathway")
}

#' This function lists all pathway components of a given pathway.
#' 
#' This function returns a (unique) data.frame listing all component IDs, names and classes of the supplied pathway.
#' 
#' @param biopax A biopax model
#' @param id string. A pathway ID
#' @param includeSubPathways logical. If TRUE the returned list will include subpathways and pathwaysteps as well.
#' @return data.frame
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  listPathwayComponents(biopax, id="pid_p_100002_wntpathway")
listPathwayComponents <- function(biopax, id, includeSubPathways=TRUE) {
	#support for bp level 3
	pwcompname = tolower("PATHWAY-COMPONENTS")
	subpathwayproperties = c("STEP-INTERACTIONS","NEXT-STEP")
	if(biopax$biopaxlevel == 3) {
		pwcompname = tolower("pathwayComponent")
		subpathwayproperties = c("nextStep","stepProcess","pathwayOrder")
	}
	#####CONTINUE HERE
	
	
	id = unique(striphash(id))
	#get pw component list
	#pwcomp_list = as.character(biopax$df[tolower(biopax$df$property) == pwcompname & biopax$df$id %in% id,"property_attr_value"])
	if(includeSubPathways)	{
		pwcomp_list = getReferencedIDs(biopax, id, recursive=TRUE, onlyFollowProperties=c(pwcompname, subpathwayproperties))
	} else {
		pwcomp_list = getReferencedIDs(biopax, id, recursive=TRUE, onlyFollowProperties=c(pwcompname))
	}
	#strip # from front of id
	listInstances(biopax, id=unique(striphash(pwcomp_list)))
}

#' This function lists all components of a given complex.
#' 
#' This function returns a (unique) data.frame listing all component IDs, names and classes of the supplied complex.
#' 
#' @param biopax A biopax model
#' @param id string. A complex ID
#' @return data.frame
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  listComplexComponents(biopax, id="ex_m_100650")
listComplexComponents <- function(biopax, id) {
	#support for bp level 3
	compname = "COMPONENTS"
	if(biopax$biopaxlevel == 3) {
		compname = "component"
	}
	
	id = unique(striphash(id))
	#get complex component list
	complexcomp_list = as.character(unique(unlist(getReferencedIDs(biopax, id, recursive=FALSE, onlyFollowProperties=c(compname)))))
	#strip # from front of id
	listInstances(biopax, id=unique(striphash(complexcomp_list)))
}

#' This function lists all components of a given interaction.
#' 
#' This function returns a (unique) data.frame listing IDs, names and classes of all components of the supplied interaction.
#' 
#' @param biopax A biopax model
#' @param id string. A complex ID
#' @param splitComplexes logical. If TRUE complexes are split up into their components and the added to the listing.
#' @return data.frame
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  listInteractionComponents(biopax, id="ex_i_100036_activator_1")
listInteractionComponents <- function(biopax, id, splitComplexes=TRUE) {
	#support for bp level 3
	if(biopax$biopaxlevel == 2) {
		compname = c("PARTICIPANTS","CONTROLLER","CONTROLLED","LEFT","RIGHT","COFACTOR","PHYSICAL-ENTITY")
		if(splitComplexes) compname = c(compname,"COMPONENTS")
	}
	if(biopax$biopaxlevel == 3) {
		compname = c("participant","controller","controlled","left","right","cofactor","product","template")
		if(splitComplexes) compname = c(compname,"component")
	}
	
	id = unique(striphash(id))
	#get complex component list
	comp_list = as.character(unique(unlist(getReferencedIDs(biopax, id, recursive=TRUE, onlyFollowProperties=c(compname)))))
	#strip # from front of id
	listInstances(biopax, id=unique(striphash(comp_list)))
}

#' This function generates the gene set of a pathway.
#'  
#' This function generates a gene set of all physicalEntity's of a pathway. First all interactions of the pathway are retrieved and all components of these interactions are then listed.  
#'  
#' @param biopax A biopax model
#' @param pwid string
#' @return Returns the gene set of the supplied pathway. Returns NULL if the pathway has no components.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  pwid1 = "pid_p_100002_wntpathway"
#'  pathway2Geneset(biopax, pwid=pwid1)
pathway2Geneset <- function(biopax, pwid) {
	
	pwComponents = listPathwayComponents(biopax, id=pwid)$id
	interactionComponents = NULL
	if(length(pwComponents)>0) {
		for(p in pwComponents) {
			interactionComponents = c(interactionComponents, listInteractionComponents(biopax,id=p)$id)
		}
		return(listInstances(biopax, id=interactionComponents))
	} else {
		return(NULL)
	}
	
}



#' This functions splits up a complex into its components.
#' 
#' This function looks up the supplied Complex ID and returns the names of all its components. 
#' 
#' @param biopax A biopax model 
#' @param complexid string ID of an complex
#' @param recursive logical
#' @return Returns a character vector with the names of all subcomponents.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  selectInstances(biopax, id="ex_m_100650")
#'  listInstances(biopax, id="ex_m_100650")
#'  listComplexComponents(biopax, id="ex_m_100650")
#'  splitComplex(biopax, complexid="ex_m_100650")
splitComplex <- function(biopax, complexid, recursive=TRUE) {
	#support for bp level 3
	compname = c("COMPONENTS","PHYSICAL-ENTITY")
	if(biopax$biopaxlevel == 3) {
		compname = c("component")
	}
	
	# complexes can contain entries via "bp:COMPONENTS" -> physicalentityparticipant "bp:PHYSICAL-ENTITY" -> physicalentity
	referenced = selectInstances(biopax, id=getReferencedIDs(biopax, complexid, recursive=recursive, onlyFollowProperties=compname))
	
	sel = tolower(referenced$class) %in% tolower(c("dna","rna","protein","smallMolecule"))
	
	#unique(referenced[tolower(referenced$property)==tolower("NAME") & sel,c("id","property_value")])
	
	listInstances(biopax,id=as.character(unique(referenced[sel,"id"])))
}


#' This function returns a vector of ids of all instances referenced by the specified instance.
#' 
#' This function takes an id and a biopax model as input. The id of every instance that is referenced is returned.
#' If recursive == TRUE this function recurses through all referenced IDs of the referenced instances and so on.
#' "onlyFollowProperties" limits the recursivness to only certain properties, for example follow only complexes or physicalEntities.
#' 
#' @param biopax A biopax model 
#' @param id string. ID of the instance
#' @param recursive logical
#' @param onlyFollowProperties character vector
#' @return Returns a character vector of IDs referenced by the supplied id in the supplied biopax model.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  listComplexComponents(biopax, id="ex_m_100650")
#'  getReferencedIDs(biopax, id="ex_m_100650", recursive=FALSE)
#'  getReferencedIDs(biopax, id="ex_m_100650", recursive=TRUE)
getReferencedIDs <- function(biopax, id, recursive=TRUE, onlyFollowProperties=c()) {
	id = unique(striphash(id))
	referencedIDs = list()
	
	##speed up selecting by doing it once:
	isrdfresource = biopax$df$property_attr == "rdf:resource"
	
	#every ref in instances of id
	if(length(onlyFollowProperties) > 0) {
		propertysel = tolower(biopax$df$property) %in% tolower(onlyFollowProperties)
		newIDs = biopax$df[biopax$df$id %in% id & isrdfresource & propertysel,"property_attr_value"]	
	} else {
		newIDs = biopax$df[biopax$df$id %in% id & isrdfresource ,"property_attr_value"]	
	}
	newIDs = unique(striphash(newIDs))
	newIDs = newIDs[!(newIDs %in% id)]
	referencedIDs = c(referencedIDs,newIDs)
	
	if(recursive) {
		while(length(newIDs)>0) {
			if(length(onlyFollowProperties) > 0) {
				newIDs = biopax$df[biopax$df$id %in% newIDs & isrdfresource & propertysel,"property_attr_value"]	
			} else {
				newIDs = biopax$df[biopax$df$id %in% newIDs & isrdfresource ,"property_attr_value"]	
			}
			newIDs = unique(striphash(newIDs))
			newIDs = newIDs[!(newIDs %in% c(referencedIDs,id))]
			referencedIDs = c(referencedIDs,newIDs)
		}
	}
	unlist(referencedIDs)
}

#' This function returns a vector of ids of all instances that reference the supplied id.
#' 
#' This function takes an id and a biopax model as input. The id of every instance that references the supplied id is returned.
#' If recursive == TRUE this function recurses through all referencing IDs of the referencing instances and so on.
#' "onlyFollowProperties" limits the recursivness to only certain properties, for example follow only complexes or physicalEntities.
#' 
#' @param biopax A biopax model 
#' @param id string. ID of the instance
#' @param recursive logical
#' @param onlyFollowProperties character vector
#' @return Returns a character vector of IDs referencing the supplied id in the supplied biopax model.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  listComplexComponents(biopax, id="ex_m_100650")
#'  getReferencingIDs(biopax, id="ex_m_100650", recursive=FALSE)
#'  getReferencingIDs(biopax, id="ex_m_100650", recursive=TRUE)
getReferencingIDs <- function(biopax, id, recursive=TRUE, onlyFollowProperties=c()) {
	id = unique(id)
	id = addhash(id)
	referencingIDs = list()
	
	##speed up selecting by doing it once:
	isrdfresource = biopax$df$property_attr == "rdf:resource"

	if(length(onlyFollowProperties) > 0) {
		propertysel = tolower(biopax$df$property) %in% tolower(onlyFollowProperties)
		newIDs = biopax$df[biopax$df$property_attr_value %in% id & isrdfresource & propertysel,"id"]	
	} else {
		newIDs = biopax$df[biopax$df$property_attr_value %in% id & isrdfresource ,"id"]	
	}
	newIDs = unique(addhash(newIDs))
	newIDs = newIDs[!(newIDs %in% id)]
	if(length(newIDs)>0) referencingIDs = c(referencingIDs, newIDs)
	
	if(recursive) {
		while(length(newIDs)>0) {
			if(length(onlyFollowProperties) > 0) {
				newIDs = biopax$df[biopax$df$property_attr_value %in% newIDs & isrdfresource & propertysel,"id"]	
			} else {
				newIDs = biopax$df[biopax$df$property_attr_value %in% newIDs & isrdfresource ,"id"]	
			}
			newIDs = unique(addhash(newIDs))
			newIDs = newIDs[!(newIDs %in% c(referencingIDs,id))]
			if(length(newIDs)>0) referencingIDs = c(referencingIDs, newIDs)
		}
	}
	unlist(unique(striphash(referencingIDs)))
}

#' This function returns the class name of the instance.
#' 
#' This function returns the class name of the instance.
#' 
#' @param biopax A biopax model
#' @param id string
#' @return Returns the class name of the biopax instance.
#' @author fkramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  getInstanceClass(biopax, id="ex_m_100650")
getInstanceClass <- function(biopax, id) {
	id = striphash(id)
	as.character(biopax$df[biopax$df$id == id,"class"][1])
}

#' This function returns all properties of the specified type for an instance.
#' 
#' This function returns all properties of the specified type for an instance. By default this function returns the NAME property of an instance.
#' 
#' @param biopax A biopax model
#' @param id string
#' @param property string. Attention: All properties in Biopax Level 2 are all upper case. 
#' @return Returns a character vector with all properties of the selected type for this instance. Returns NULL if no property of this type is found.
#' @author fkramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  getInstanceProperty(biopax, id="ex_m_100650", property="NAME")
#'  getInstanceProperty(biopax, id="ex_m_100650", property="ORGANISM")
#'  getInstanceProperty(biopax, id="ex_m_100650", property="COMPONENTS")
getInstanceProperty <- function(biopax, id, property="NAME") {
	id = striphash(id)
	
	### speed up and quick fix for biopax level 3 naming:
	if(tolower(property) == "name") {
		#names = unfactorize(selectInstances(biopax, id=id, property=))
		#names = biopax$df[biopax$df$id==id & tolower(biopax$df$property) %in% c("name","displayname","standardname"),]
		names = biopax$df[biopax$df$id==id,]
		if(dim(names)[1] == 0) return(NULL)
		names = names[tolower(names$property) %in% c("name","displayname","standardname"),]
		if(dim(names)[1] == 0) return(NULL)
		if(any(grepl("displayName", names$property))) return(as.character(names[names$property == "displayName", "property_value"]))
		if(any(grepl("standardName", names$property))) return(as.character(names[names$property == "standardName", "property_value"]))
		if(any(grepl("name", names$property, ignore.case=T))) return(as.character(names[tolower(names$property) == "name", "property_value"]))
		return(NULL)
	}
	
	#get class of the instance
	class = getInstanceClass(biopax,id)
	#get the properties this instance has
	classproperties = getClassProperties(class)
	#if its a reference return the property_attr_value (its the referenced ids) or return property_value otherwise
	property_type = unlist(classproperties[classproperties$property == property,]$property_type)[1]
	
	#value
	if(length(property_type)>0) {
		if(grepl("string",property_type) || grepl("double",property_type) || grepl("float",property_type) || grepl("integer",property_type)) {
			as.character(biopax$df[biopax$df$id == id & tolower(biopax$df$property) == tolower(property),"property_value"])
		} else {
			as.character(biopax$df[biopax$df$id == id & tolower(biopax$df$property) == tolower(property),"property_attr_value"])
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
#' @return Returns ids of physicalEntity corresponding to the specified  physicalEntityParticipantIDs
#' @author fkramer
internal_resolvePhysicalEntityParticipant <- function(biopax, physicalEntityId) {
	unlist(getReferencedIDs(biopax, id=physicalEntityId, recursive=FALSE, onlyFollowProperties=c("PHYSICAL-ENTITY")))
}

#' This function returns the neighborhood of a physicalEntity
#' 
#' This function searches the supplied biopax for interactions that are connected to the molecule or within 'depth' number of steps from it.
#' 
#' @param biopax A biopax model
#' @param id string. ID of a physicalEntity (dna, rna, protein, complex, smallMolecule)
#' @param depth integer. Search depth, this specifies how far out from the specified molecule the neighborhood should be streched.
#' @param onlyInPathways character vector of pathway IDs. Search only in these pathways for neighbors. 
#' @return Returns ids of interactions within 'depth' number of steps of the specified physicalEntity 
#' @author fkramer
#' @export
getNeighborhood <- function(biopax, id, depth=1, onlyInPathways=c()) {
	### TODO doesnt work yet. fix getreferencinginstances, add together interactionlist + moleculelist
	
	#TODO ####################################continue here! how to traverse backwards through biopax networks??
	
	## get outgoing edges, consider complexes?
	## get instances referencing the PEP(s)/complexe(s) in CONTROLLER/PARTICIPANTS
	
	## get incomming edges. 
	## get instances with id occuring as PARTICIPANTS/LEFT/RIGHT/... in CONTROLLEDs ?!
	
	interactionlist = c()
	moleculelist = id
	for(i in 1:depth) {
		# add all interaction that includes at least one of the molecules to the interactionlist
		#biopax$df[biopax$df$id == id & biopax$df$property == property,"property_attr_value"]
		#getReferencingInstances(biopax, id, recursive=TRUE, onlyFollowProperties=c("COMPONENTS","PHYSICAL-ENTITY"))
		selectInstances(biopax, id=getReferencingIDs(biopax, id=id, recursive=TRUE, onlyFollowProperties=c("COMPONENTS","PHYSICAL-ENTITY","LEFT","RIGHT","CONTROLLER","CONTROLLED")))
	}
}

#' This function returns the annotations of the supplied instances.
#' 
#' This function returns the annotations of the supplied IDs in a data.frame.
#' 
#' @param biopax A biopax model
#' @param id vector of strings. IDs of instances to get annotations
#' @param splitComplexes logical. If TRUE complexes are split up into their components and the annotation of the components is added.
#' @param followPhysicalEntityParticipants logical. If TRUE physicalEntityParticipants are resolved to their corresponding physicalEntities and their annotation is added. 
#' @return Returns data.frame with annotations 
#' @author fkramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#' # no annotations for exactly the complex
#' getXrefAnnotations(biopax, id="ex_m_100650")
#' # split up the complex and get annotations for all the molecules involved
#' getXrefAnnotations(biopax, id="ex_m_100650", splitComplexes=TRUE)
getXrefAnnotations <- function(biopax, id, splitComplexes=FALSE, followPhysicalEntityParticipants=TRUE) {
	id = c(striphash(id))
	#annotations = data.frame(type=NA, id=NA, name=NA, annotation_type=NA, annotation_id=NA, annotation=NA, stringsAsFactors=FALSE)
	annotations = matrix(nrow=0, ncol=6)
	colnames(annotations) = c("type", "id","name","annotation_type", "annotation_id", "annotation")
	
	for(i in 1:length(id)) {
		# if its a complex AND we're supposed to split it:
		if(tolower(getInstanceClass(biopax,id[i])) == "complex" & splitComplexes) {
			referenced = selectInstances(biopax, id=getReferencedIDs(biopax, id[i], recursive=TRUE, onlyFollowProperties=c("COMPONENTS","PHYSICAL-ENTITY")))
			sel = tolower(referenced$class) %in% tolower(c("dna","rna","protein","smallMolecule"))
			referenced = as.character(unique(referenced[sel & !(referenced$id %in% id),"id"]))
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

