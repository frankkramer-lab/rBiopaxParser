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
#' Returns all instances that conform to the selection criteria. This function returns a subset of the internal data.table of the biopax object.
#' Selection criteria are wether instances belong to a certain class or have the specified id, property or name. Setting a criteria to NULL ignores this criteria.
#' If returnValues is set to FALSE only the selector (a logical vector with length of the internal data.table) is returned, otherwise the selected data is returned.
#' If includeSubClasses is set to TRUE the class criteria is broadened to include all classes that inherit from the given class, e.g. if class="control" and includeSubClasses=TRUE the function will select catalyses and modulations too, since they are a subclass of class control. 
#' If includeReferencedInstances is set to TRUE all instances that are being referenced by the selected instances are being selected too. The parameter works recursively, this means for example that a selected pathway and all it's interactions, complexes, molecules and annotations are returned if this parameter is set to true. This parameter is especially helpful if you want to migrate or merge knowledge from different data bases.
#' 
#' @param biopax A biopax model or a compatible internal data.table
#' @param id string. ID of the instances to select
#' @param class string. Class of the instances to select
#' @param property string. Return only this property of the instances
#' @param name string. Name of the instances to select
#' @param returnValues logical. If returnValues is set to FALSE only the selector (a logical vector with length of the internal data.table) is returned, otherwise the selected data is returned
#' @param includeSubClasses logical. If includeSubClasses is set to TRUE the class criteria is broadened to include all classes that inherit from the given class
#' @param includeReferencedInstances logical. If includeReferencedInstances is set to TRUE all instances that are being referenced by the selected instances are being selected too
#' @param returnCopy logical. Defaults to TRUE. If TRUE a copy of the internal data.table is returned. If FALSE data is returned by reference. Set to FALSE to increase speed when only ever reading data. Make sure you understand the implications of using this! See vignette of data.table package.
#' @param biopaxlevel integer. Set the biopax level here if you supply a data.table directly.
#' @return Returns a data.table containing all instances conforming to the given selection criteria if returnValues=TRUE, only the selector for the internal data.table otherwise.
#' @author Frank Kramer
#' @export
#' @import data.table
#' @examples
#'  # load data
#'  data(biopax2example)
#'  # select the subset of the internal data.table that belongs to class "protein"
#'  selectInstances(biopax, class="protein")
#'  # select the subset of the internal data.table that belongs to class "interaction"
#'  selectInstances(biopax, class="interaction")
#'  # select the subset of the internal data.table that belongs to class "interaction" or any of its sub classes, like control, catalysis etc.
#'  selectInstances(biopax, class="interaction", includeSubClasses=TRUE)
#'  # select the subset of the internal data.table that belongs to class "pathway" AND is a "NAME" property
#'  selectInstances(biopax, class="pathway", property="NAME")
selectInstances <- function (biopax, id=NULL, class=NULL, property=NULL, name=NULL, returnValues=TRUE, includeSubClasses=FALSE, includeReferencedInstances=FALSE, returnCopy=TRUE, biopaxlevel=NULL) {
	var_id=id
	rm(id)
	var_class = class
	rm(class)
	
	if("biopax" %in% class(biopax)) {
		df = biopax$dt
		biopaxlevel = biopax$biopaxlevel
	} else if("biopax_df" %in% class(biopax)) {
		df = biopax
		if(is.null(biopaxlevel))  biopaxlevel=3
	}  else {
		stop("selectInstances: parameter biopax is neither biopax object nor compatible biopax data.table")
	}
	
	sel = rep(TRUE,length.out=dim(df)[1])
	
	if(!is.null(var_id)) {
		var_id = unique(striphash(var_id))
		sel = sel & (df$id %chin% var_id)
	}
	
	if(!is.null(var_class)) {
		if(includeSubClasses) {
			var_class = unique(c(var_class,getSubClasses(var_class, biopaxlevel)))
		}
		sel = sel & (tolower(df$class) %chin% tolower(stripns(var_class)))
	}
	
	if(!is.null(property)) {
		sel = sel & (tolower(df$property) %chin% tolower(property))
	}
	
	if(!is.null(name)) {
		ids = as.character(df[property_value %chin% name][tolower(property)=="name"]$id)
		ids = unique(ids)
		sel = sel & (df$id %chin% ids)
	}
	
	if(includeReferencedInstances) {
		#include all referenced instances. this is the only place we do a logical OR.
		ids = as.character(df[sel]$id)
		ids = unique(ids)
		ids = getReferencedIDs(df, ids)
		sel = sel | (df$id %chin% ids)
	}	
	
	if(!returnValues) return(sel)
	if(returnCopy) return(copy(df[sel]))
	df[sel]
}


#' Lists all instances that conform to the selection criteria.
#' 
#' Lists all instances that conform to the selection criteria. In contrast to selectInstances this function returns an easier to read list.
#' This function returns an ordered data.table of class, id and name of the instances.
#' Selection criteria are wether instances belong to a certain class or have the specified id or name. Setting a criteria to NULL ignores this criteria.
#' If includeSubClasses is set to TRUE the class criteria is broadened to include all classes that inherit from the given class, e.g. if class="control" and includeSubClasses=TRUE the function will select catalyses and modulations too, since they are a subclass of class control. 
#' 
#' @param biopax A biopax model
#' @param id string. ID of the instances to select
#' @param class string. Class of the instances to select
#' @param name string. Name of the instances to select
#' @param includeSubClasses logical. If includeSubClasses is set to TRUE the class criteria is broadened to include all classes that inherit from the given class
#' @param returnIDonly logical. If TRUE only IDs of the components are returned. This saves time for looking up names for every single ID.
#' @param biopaxlevel integer. Set the biopax level here if you supply a data.table directly.
#' @return Returns a data.frame containing all instances conforming to the given selection criteria. If returnIDonly=TRUE, only the selector for the internal data.table otherwise.
#' @author Frank Kramer
#' @import data.table
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
listInstances <- function (biopax, id=NULL, class=NULL, name=NULL, includeSubClasses=FALSE, returnIDonly=FALSE, biopaxlevel=NULL) {
	var_id=id
	rm(id)
	var_class = class
	rm(class)
	
	if("biopax" %in% class(biopax)) {
		bpsel = biopax$dt
		biopaxlevel = biopax$biopaxlevel
	} else if("biopax_df" %in% class(biopax)) {
		bpsel = biopax
		if(is.null(biopaxlevel)) biopaxlevel=3
	}  else {
		stop("listInstances: parameter biopax is neither biopax object nor compatible biopax data.table")
	}
	
	if(!is.null(var_id)) {
		var_id = unique(striphash(var_id))
		bpsel = bpsel[id %chin% var_id]
	}
	
	if(!is.null(var_class)) {
		if(includeSubClasses) {
			var_class = unique(c(var_class,getSubClasses(var_class, biopaxlevel)))
		}
		bpsel = bpsel[tolower(class) %chin% tolower(stripns(var_class))]
	}
	
	if(!is.null(name)) {
		ids = as.character(bpsel[property_value %chin% name][tolower(property) %chin% c("name","displayname","standardname")]$id)
		ids = unique(ids)
		bpsel = bpsel[id %chin% ids]
	}
	
	# exceptions: empty selection & return only ids
	if(dim(bpsel)[1]==0) return(NULL)
	if(returnIDonly) return(unique(bpsel$id))
	
	props = tolower(bpsel$property)
	names = c(which(props=="displayname", arr.ind=T),which(props=="standardname", arr.ind=T),which(props=="name", arr.ind=T))
	names = bpsel[names, list(id,property_value)]
	names = names[!duplicated(names,by="id")]
	extraids =  bpsel[!(id %chin% names$id)]$id
	if(length(extraids) >0)	names = rbindlist(list(names, data.table(id=extraids, property_value="")))
	setnames(names,c("id","name"))
	ret = as.data.frame(names)
	rm(bpsel)
	return(ret)
	
}


#' This function returns a list of all pathway ids.
#' 
#' This function returns a vector of all pathway ids.
#' 
#' @param biopax A biopax model 
#' @return Returns a character vector containing the names of all pathways.
#' @author Frank Kramer
#' @export
#' @import data.table
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
#' @param returnIDonly logical. If TRUE only IDs of the components are returned. This saves tiem for looking up names for every single ID.
#' @return data.frame
#' @author Frank Kramer
#' @export
#' @import data.table
#' @examples
#'  # load data
#'  data(biopax2example)
#'  listPathwayComponents(biopax, id="pid_p_100002_wntpathway")
listPathwayComponents <- function(biopax, id, includeSubPathways=TRUE, returnIDonly=FALSE) {
	#support for bp level 2 and 3
	if(biopax$biopaxlevel == 2) {
		pwcompname = "pathway-components"
		subpathwayproperties = c("step-interactions","next-step")
	}
	if(biopax$biopaxlevel == 3) {
		pwcompname = "pathwaycomponent"
		subpathwayproperties = c("nextStep","stepProcess","pathwayOrder")
	}
	
	id = unique(striphash(id))
	id = id[!is.na(id) & !is.null(id) & nchar(id) > 0 ]
	#get pw component list
	#pwcomp_list = as.character(biopax$dt[tolower(biopax$dt$property) == pwcompname & biopax$dt$id %in% id,"property_attr_value"])
	if(includeSubPathways)	{
		pwcomp_list = getReferencedIDs(biopax, id, recursive=TRUE, onlyFollowProperties=c(pwcompname, subpathwayproperties))
	} else {
		pwcomp_list = getReferencedIDs(biopax, id, recursive=TRUE, onlyFollowProperties=c(pwcompname))
	}
	if(is.null(pwcomp_list)) return(NULL)
	pwcomp_list = pwcomp_list[!is.na(pwcomp_list) & !is.null(pwcomp_list) & nchar(pwcomp_list) > 0 ]
	if(returnIDonly) return(unique(striphash(pwcomp_list)))
	#strip # from front of id
	listInstances(biopax, id=unique(striphash(pwcomp_list)))
}

#' This function lists all components of a given complex.
#' 
#' This function returns a (unique) data.frame listing all component IDs, names and classes of the supplied complex.
#' 
#' @param biopax A biopax model
#' @param id string. A complex ID
#' @param returnIDonly logical. If TRUE only IDs of the components are returned. This saves tiem for looking up names for every single ID.
#' @return data.frame
#' @author Frank Kramer
#' @export
#' @import data.table
#' @examples
#'  # load data
#'  data(biopax2example)
#'  listComplexComponents(biopax, id="ex_m_100650")
listComplexComponents <- function(biopax, id, returnIDonly=FALSE) {
	#support for bp level 3
	compname = "COMPONENTS"
	if(biopax$biopaxlevel == 3) {
		compname = "component"
	}
	
	id = unique(striphash(id))
	#get complex component list
	complexcomp_list = as.character(unique(getReferencedIDs(biopax, id, recursive=FALSE, onlyFollowProperties=c(compname))))
	if(is.null(complexcomp_list)) return(NULL)
	complexcomp_list = complexcomp_list[!is.na(complexcomp_list) & !is.null(complexcomp_list) & nchar(complexcomp_list) > 0 ]
	#strip # from front of id
	if(returnIDonly) return(unique(striphash(complexcomp_list)))
	listInstances(biopax, id=unique(striphash(complexcomp_list)))
}

#' This function lists all components of a given interaction.
#' 
#' This function returns a (unique) data.frame listing IDs, names and classes of all components of the supplied interaction.
#' 
#' @param biopax A biopax model
#' @param id string. A complex ID
#' @param splitComplexes logical. If TRUE complexes are split up into their components and the added to the listing.
#' @param returnIDonly logical. If TRUE only IDs of the components are returned. This saves tiem for looking up names for every single ID.
#' @return data.frame
#' @author Frank Kramer
#' @export
#' @import data.table
#' @examples
#'  # load data
#'  data(biopax2example)
#'  listInteractionComponents(biopax, id="ex_i_100036_activator_1")
listInteractionComponents <- function(biopax, id, splitComplexes=TRUE, returnIDonly=FALSE) {
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
	if(is.null(comp_list)) return(NULL)
	comp_list = comp_list[!is.na(comp_list) & !is.null(comp_list) & nchar(comp_list) > 0 ]
	#strip # from front of id
	if(returnIDonly) return(unique(striphash(comp_list)))
	listInstances(biopax, id=unique(striphash(comp_list)))
}

#' This function generates the gene set of a pathway.
#'  
#' This function generates a gene set of all physicalEntity's of a pathway. First all interactions of the pathway are retrieved and all components of these interactions are then listed.  
#'  
#' @param biopax A biopax model
#' @param pwid string
#' @param returnIDonly logical. If TRUE only IDs of the components are returned. This saves tiem for looking up names for every single ID.
#' @return Returns the gene set of the supplied pathway. Returns NULL if the pathway has no components.
#' @author Frank Kramer
#' @export
#' @import data.table
#' @examples
#'  # load data
#'  data(biopax2example)
#'  pwid1 = "pid_p_100002_wntpathway"
#'  pathway2Geneset(biopax, pwid=pwid1)
pathway2Geneset <- function(biopax, pwid, returnIDonly=FALSE) {

	pwComponents = listPathwayComponents(biopax, id=pwid, returnIDonly=TRUE)
	interactionComponents = NULL
	if(length(pwComponents)>0) {
		for(p in pwComponents) {
			interactionComponents = c(interactionComponents, listInteractionComponents(biopax,id=p, returnIDonly=T))
		}
		interactionComponents = interactionComponents[!is.na(interactionComponents) & !is.null(interactionComponents) & nchar(interactionComponents) > 0 ]
		if(is.null(interactionComponents)) return(NULL)
		if(returnIDonly) return(unique(interactionComponents))
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
#' @param returnIDonly logical. If TRUE only IDs of the components are returned. This saves tiem for looking up names for every single ID.
#' @param biopaxlevel integer. Set the biopax level here if you supply a data.table directly.
#' @return Returns a character vector with the names of all subcomponents.
#' @author Frank Kramer
#' @export
#' @import data.table
#' @examples
#'  # load data
#'  data(biopax2example)
#'  selectInstances(biopax, id="ex_m_100650")
#'  listInstances(biopax, id="ex_m_100650")
#'  listComplexComponents(biopax, id="ex_m_100650")
#'  splitComplex(biopax, complexid="ex_m_100650")
splitComplex <- function(biopax, complexid, recursive=TRUE, returnIDonly=FALSE, biopaxlevel=3) {
	
	if("biopax" %in% class(biopax)) {
		df = biopax$dt
		biopaxlevel = biopax$biopaxlevel
	} else if("biopax_df" %in% class(biopax)) {
		df = biopax
	}  else {
		stop("splitComplex: parameter biopax is neither biopax object nor compatible biopax data.table")
	}
	
	#support for bp level 3
	compname = c("COMPONENTS","PHYSICAL-ENTITY")
	if(biopaxlevel == 3) {
		compname = c("component")
	}
	
	# complexes can contain entries via "bp:COMPONENTS" -> physicalentityparticipant "bp:PHYSICAL-ENTITY" -> physicalentity
	ref = getReferencedIDs(df, complexid, recursive=recursive, onlyFollowProperties=compname)
	if(is.null(ref)) return(NULL)
	referenced = selectInstances(df, id=ref, returnCopy=FALSE)
	
	referenced = unique(as.character(referenced[tolower(class) %chin% c("dna","rna","protein","smallmolecule")]$id))
	if(length(referenced)==0) return(NULL)
	if(returnIDonly) return(striphash(referenced))
	
	listInstances(df,id=referenced)
}


#' This function returns a vector of ids of all instances referenced by the specified instance.
#' 
#' This function takes an id and a biopax model as input. The id of every instance that is referenced is returned.
#' If recursive == TRUE this function recurses through all referenced IDs of the referenced instances and so on.
#' "onlyFollowProperties" limits the recursivness to only certain properties, for example follow only complexes or physicalEntities.
#' 
#' @param biopax A biopax model OR a compatible data.table
#' @param id string. ID of the instance
#' @param recursive logical
#' @param onlyFollowProperties character vector
#' @return Returns a character vector of IDs referenced by the supplied id in the supplied biopax model.
#' @author Frank Kramer
#' @export
#' @import data.table
#' @examples
#'  # load data
#'  data(biopax2example)
#'  listComplexComponents(biopax, id="ex_m_100650")
#'  getReferencedIDs(biopax, id="ex_m_100650", recursive=FALSE)
#'  getReferencedIDs(biopax, id="ex_m_100650", recursive=TRUE)
getReferencedIDs <- function(biopax, id, recursive=TRUE, onlyFollowProperties=c()) {
	var_id = unique(striphash(id))
	rm(id)
	referencedIDs = vector()
	
	if("biopax" %in% class(biopax)) {
		bpsel = biopax$dt[property_attr == "rdf:resource",]
	} else if("biopax_df" %in% class(biopax)) {
		bpsel = biopax[property_attr == "rdf:resource",]
	}  else {
		stop("getReferencedIDs: parameter biopax is neither biopax object nor compatible biopax data.table")
	}
	
	#every ref in instances of id
	if(length(onlyFollowProperties) > 0) {
		propertysel = tolower(bpsel$property) %chin% tolower(onlyFollowProperties)
		newIDs = bpsel[propertysel][id %chin% var_id]$property_attr_value	
	} else {
		newIDs =  bpsel[id %chin% var_id]$property_attr_value
	}
	
	if(length(newIDs)==0) return(NULL)
	newIDs = unique(striphash(newIDs))
	newIDs = newIDs[!(newIDs %chin% var_id)]
	referencedIDs = c(referencedIDs,newIDs)
	
	if(recursive) {
		while(length(newIDs)>0) {
			if(length(onlyFollowProperties) > 0) {
				newIDs = bpsel[propertysel][id %chin% newIDs]$property_attr_value
			} else {
				newIDs = bpsel[id %chin% newIDs]$property_attr_value
			}
			newIDs = unique(striphash(newIDs))
			newIDs = newIDs[!(newIDs %chin% c(referencedIDs,var_id))]
			referencedIDs = c(referencedIDs,newIDs)
		}
	}
	
	ret = unique(striphash(referencedIDs))
	ret = ret[ret!=""]
	if(length(ret)==0) return(NULL)
	return(ret)

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
#' @import data.table
#' @examples
#'  # load data
#'  data(biopax2example)
#'  listComplexComponents(biopax, id="ex_m_100650")
#'  getReferencingIDs(biopax, id="ex_m_100650", recursive=FALSE)
#'  getReferencingIDs(biopax, id="ex_m_100650", recursive=TRUE)
getReferencingIDs <- function(biopax, id, recursive=TRUE, onlyFollowProperties=c()) {
	var_id = unique(id)
	rm(id)
	var_id = addhash(var_id)
	referencingIDs = vector()
	
	if("biopax" %in% class(biopax)) {
		bpsel = biopax$dt[property_attr == "rdf:resource",]
	} else if("biopax_df" %in% class(biopax)) {
		bpsel = biopax[property_attr == "rdf:resource",]
	}  else {
		stop("getReferencedIDs: parameter biopax is neither biopax object nor compatible biopax data.table")
	}
	

	if(length(onlyFollowProperties) > 0) {
		propertysel = tolower(bpsel$property) %chin% tolower(onlyFollowProperties)
		newIDs = bpsel[propertysel][property_attr_value %chin% var_id]$id
	} else {
		newIDs = bpsel[property_attr_value %chin% var_id]$id
	}
	
	if(length(newIDs)==0) return(NULL)
	newIDs = unique(addhash(newIDs))
	newIDs = newIDs[!(newIDs %chin% var_id)]
	if(length(newIDs)>0) referencingIDs = c(referencingIDs, newIDs)
	
	if(recursive) {
		while(length(newIDs)>0) {
			if(length(onlyFollowProperties) > 0) {
				newIDs = bpsel[propertysel][property_attr_value %chin% newIDs]$id
			} else {
				newIDs = bpsel[property_attr_value %chin% newIDs]$id
			}
			newIDs = unique(addhash(newIDs))
			newIDs = newIDs[!(newIDs %chin% c(referencingIDs,var_id))]
			if(length(newIDs)>0) referencingIDs = c(referencingIDs, newIDs)
		}
	}

	ret = unique(striphash(referencingIDs))
	ret = ret[ret!=""]
	if(length(ret)==0) return(NULL)
	return(ret)
	
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
#' @import data.table
#' @examples
#'  # load data
#'  data(biopax2example)
#'  getInstanceClass(biopax, id="ex_m_100650")
getInstanceClass <- function(biopax, id) {
	var_id = striphash(id)
	rm(id)
	
	if("biopax" %in% class(biopax)) {
		df = biopax$dt
	} else if("biopax_df" %in% class(biopax)) {
		df = biopax
	}  else {
		stop("getInstanceClass: parameter biopax is neither biopax object nor compatible biopax data.table")
	}
	
	as.character(df[id %chin% var_id]$class[1])
}

#' This function returns all properties of the specified type for an instance.
#' 
#' This function returns all properties of the specified type for an instance. By default this function returns the NAME property of an instance.
#' 
#' @param biopax A biopax model
#' @param id string
#' @param property string.
#' @param includeAllNames logical. Biopax Level 3 brought 2 new name properties: displayName and standardName. Per default this return all names of an instance. Disable if you only want the NAME property.  
#' @return Returns a character vector with all properties of the selected type for this instance. Returns NULL if no property data is found.
#' @author fkramer
#' @export
#' @import data.table
#' @examples
#'  # load data
#'  data(biopax2example)
#'  getInstanceProperty(biopax, id="ex_m_100650", property="NAME")
#'  getInstanceProperty(biopax, id="ex_m_100650", property="ORGANISM")
#'  getInstanceProperty(biopax, id="ex_m_100650", property="COMPONENTS")
getInstanceProperty <- function(biopax, id, property="NAME", includeAllNames=TRUE) {
	if(is.null(id) | is.na(id)) return(NULL)
	var_id = striphash(id)
	rm(id)
	var_property = tolower(property)
	rm(property)
	
	if("biopax" %in% class(biopax)) {
		df = biopax$dt
	} else if("biopax_df" %in% class(biopax)) {
		df = biopax
	}  else {
		stop("getInstanceProperty: parameter biopax is neither biopax object nor compatible biopax data.table")
	}
	
	### speed up and quick fix for biopax level 3 naming:
	if(var_property == "name") {
		names = df[id==var_id,]
		if(dim(names)[1] == 0) return(NULL)
		if(includeAllNames) {
			names = names[tolower(property) %chin% c("name","displayname","standardname"),] 
		} else {
			names = names[tolower(property) == "name",]
		}
		if(dim(names)[1] == 0) return(NULL)
		return(as.character(names$property_value))
	}
	
	data = selectInstances(df, id=var_id, returnCopy = FALSE)
	data = data[tolower(property) %chin% var_property,]
	
	#value
	if(dim(data)[1]>0) {
		if(grepl("string",data$property_attr_value) || grepl("double",data$property_attr_value) || grepl("float",data$property_attr_value) || grepl("integer",data$property_attr_value)) {
			return(as.character(data$property_value))
		} else {
			return(as.character(data$property_attr_value))
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
#' @import data.table
internal_resolvePhysicalEntityParticipant <- function(biopax, physicalEntityId) {
	getReferencedIDs(biopax, id=physicalEntityId, recursive=FALSE, onlyFollowProperties=c("PHYSICAL-ENTITY"))
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
#' @import data.table
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
		#biopax$dt[biopax$dt$id == id & biopax$dt$property == property,"property_attr_value"]
		#getReferencingInstances(biopax, id, recursive=TRUE, onlyFollowProperties=c("COMPONENTS","PHYSICAL-ENTITY"))
		selectInstances(biopax, id=getReferencingIDs(biopax, id=id, recursive=TRUE, onlyFollowProperties=c("COMPONENTS","PHYSICAL-ENTITY","LEFT","RIGHT","CONTROLLER","CONTROLLED")))
	}
}

#' This function returns the annotations of the supplied instances.
#' 
#' This function returns the annotations of the supplied IDs in a data.table.
#' 
#' @param biopax A biopax model
#' @param id vector of strings. IDs of instances to get annotations
#' @param splitComplexes logical. If TRUE complexes are split up into their components and the annotation of the components is added.
#' @param followPhysicalEntityParticipants logical. If TRUE physicalEntityParticipants are resolved to their corresponding physicalEntities and their annotation is added. 
#' @return Returns data.table with annotations 
#' @author fkramer
#' @export
#' @import data.table
#' @examples
#'  # load data
#'  data(biopax2example)
#' # example of annotation for a protein:
#' getXrefAnnotations(biopax, id="ex_m_100647")
#' # no annotations for exactly the complex
#' getXrefAnnotations(biopax, id="ex_m_100650")
#' # split up the complex and get annotations for all the molecules involved
#' getXrefAnnotations(biopax, id="ex_m_100650", splitComplexes=TRUE)
getXrefAnnotations <- function(biopax, id, splitComplexes=FALSE, followPhysicalEntityParticipants=TRUE) {
	if(length(id)==0) return(NULL)
	id = c(unique(striphash(id)))
	var_id=id
	rm(id)

	
	annotations = data.table(class="",id="", name="", annotation_type="", annotation_id="", annotation="")[0]
	
	bpsel = selectInstances(biopax, var_id, includeReferencedInstances = T, returnCopy = FALSE)
	
	for(i in 1:length(var_id)) {
		instanceclass = bpsel[id==var_id[i]]$class[1]
		if(is.na(instanceclass) | is.null(instanceclass)) next;
		# if its a complex AND we're supposed to split it:
		if(splitComplexes && tolower(instanceclass) == "complex") {
			#split complex
			ref = getReferencedIDs(bpsel, var_id[i], recursive=TRUE, onlyFollowProperties=c("COMPONENTS","PHYSICAL-ENTITY","component"))
			if(is.null(ref)) next;
			referenced = bpsel[id %chin% ref,]
			referenced = referenced[tolower(class) %chin% c("dna","rna","protein","smallmolecule"),]
			referenced = as.character(unique(referenced[!(id %chin% var_id)]$id))
			if(is.null(referenced) || length(referenced)[1]==0) next;
			annotations = rbindlist(list(annotations, getXrefAnnotations(biopax,referenced, splitComplexes=FALSE, followPhysicalEntityParticipants=FALSE)))
		} else if(followPhysicalEntityParticipants && instanceclass == "physicalEntityParticipant") {
			# if its a physicalEntityParticipant
			peID = internal_resolvePhysicalEntityParticipant(bpsel, var_id[i])
			if(!is.null(peID) && length(peID)>0 && !(peID %chin% var_id)) {
				annotations = rbindlist(list(annotations, getXrefAnnotations(biopax, peID, splitComplexes=splitComplexes, followPhysicalEntityParticipants=TRUE)))
			}
		} else {
			# for any other class do this
			#get instance name
			#name = getInstanceProperty(biopax, id[i], property="NAME")[1]
			name = getInstanceProperty(bpsel, var_id[i], property="NAME")[1]
			if(is.null(name)) name=""
			#xrefs = getInstanceProperty(biopax, id[i], property="XREF")
			xrefs = getInstanceProperty(bpsel, var_id[i], property="XREF")
			if(is.null(xrefs)) xrefs = NA
			
			# if its a physicalentity AND have a BP3 entityReference: add these annotations as well!
			if(tolower(instanceclass) %chin% c("dna","dnaregion","rna","rnaregion","protein","smallmolecule")) {
				sel = getReferencedIDs(bpsel, var_id[i], onlyFollowProperties=c("entityReference","memberEntityReference","memberPhysicalEntity"))
				if(!is.null(sel)) {
					for(rerId in sel) {
						xrefs = c(xrefs, getInstanceProperty(bpsel, rerId, property="XREF"))
					}
				}
			}
			
			xrefs = xrefs[!is.na(xrefs) & !is.null(xrefs) & nchar(xrefs) > 0 ]
			if(is.null(xrefs) || length(xrefs) == 0) next;
			for(xref in xrefs) {
				
				annotations = rbindlist(list(annotations, data.table(class=instanceclass, id=var_id[i], name=name, annotation_type=bpsel[id==xref,class][1], annotation_id=striphash(xref),
								annotation=paste(getInstanceProperty(bpsel, xref, property="DB"),	":", getInstanceProperty(bpsel, xref, property="ID"), sep="")
								)))
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

