###############################################################################
#
# parseBiopax.R: 	This file contains the all functions related to parsing a Biopax model.
# author: Frank Kramer <dev@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################

#' This function creates a new Biopax model from scratch
#' 
#' This function creates a new Biopax model from scratch. This is not necessary if you want to parse a BioPAX export from a file, please see: readBiopax.
#' Returns a biopax model, which is a list with named elements: 
#' \describe{
#' 		\item{df}{The data.frame representing the biopax in R}
#' 		\item{ns_rdf}{RDF Namespace}
#'		\item{ns_owl}{OWL Namespace}
#'		\item{ns_bp}{Biopax Namespace}
#' 		\item{file}{NULL}
#'	}
#'  
#' @param level integer. Specifies the BioPAX level. 
#' @return A biopax model
#' @author Frank Kramer
#' @export
#' @import data.table
#' @examples
#'  biopax = createBiopax(level=2)
#' 
createBiopax <- function(level = 3)  {
	ret = list(file = NULL)
	class(ret) <- c("biopax",class(ret))
	
	if(level==1) message("BioPAX level 1 OWL is unfortunatly not supported.")
	if(level==2) {
		ret$biopaxlevel = 2
		ret$dt = data.table(rowcount=1:5,class="",id="",property="",property_attr="",property_attr_value="",property_value="", key=c("id","class","property"))
		ret$dt = copy(ret$dt[,2:7, with=F])
		ret$dt = copy(ret$dt[0])
		class(ret$dt) = c("biopax_df",class(ret$dt))
		setkey(ret$dt, id, class, property)

		# deal with namespaces, we're interested in owl, rdf, biopax
		ret$ns_rdf = "rdf"
		ret$ns_owl = "owl"
		ret$ns_bp = "bp"
		ret$namespaces = namespaces = list(
				'rdf'="http://www.w3.org/1999/02/22-rdf-syntax-ns#",
				'bp'="http://www.biopax.org/release/biopax-level2.owl#",
				'rdfs'="http://www.w3.org/2000/01/rdf-schema#",
				'owl'="http://www.w3.org/2002/07/owl#",
				'xsd'="http://www.w3.org/2001/XMLSchema#"
		)
	} 
	if(level==3) {
		ret$biopaxlevel = 3
		ret$dt = data.table(rowcount=1:5,class="",id="",property="",property_attr="",property_attr_value="",property_value="", key=c("id","class","property"))
		ret$dt = copy(ret$dt[,2:7, with=F])
		ret$dt = copy(ret$dt[0])
		class(ret$dt) = c("biopax_df",class(ret$dt))
		setkey(ret$dt, id, class, property)
		
		# deal with namespaces, we're interested in owl, rdf, biopax
		ret$ns_rdf = "rdf"
		ret$ns_owl = "owl"
		ret$ns_bp = "bp"
		ret$namespaces = namespaces = list(
				'rdf'="http://www.w3.org/1999/02/22-rdf-syntax-ns#",
				'bp'="http://www.biopax.org/release/biopax-level3.owl#",
				'rdfs'="http://www.w3.org/2000/01/rdf-schema#",
				'owl'="http://www.w3.org/2002/07/owl#",
				'xsd'="http://www.w3.org/2001/XMLSchema#"
		)
	}
	
	ret
}


#' This function reads in a Biopax .owl file
#' 
#' This function reads in a Biopax .owl file and generates the internal data.frame format used in this package.
#' This function can take a while with really big Biopax files like NCIs Pathway Interaction Database or Reactome.
#' In almost every case this is your starting point.
#' Returns a biopax model, which is a list with named elements: 
#' \describe{
#' 		\item{df}{The data.frame representing the biopax in R}
#' 		\item{ns_rdf}{RDF Namespace}
#'		\item{ns_owl}{OWL Namespace}
#'		\item{ns_bp}{Biopax Namespace}
#' 		\item{file}{File name}
#'	}
#'  
#' @param file string. File name 
#' @param verbose logical. Output messages about how parsing is going and so on.
#' @return A biopax model
#' @author Frank Kramer
#' @export
#' @import data.table
#' @examples
#'  \dontrun{biopax = readBiopax(file="biopaxmodel.owl")}
#'  \dontrun{biopax} 
#' 
readBiopax <- function(file, verbose=TRUE)  {
	
	if(is.null(file) || !file.exists(file)) { 
		stop(paste("readBiopax: Cannot find file:",file))
		return(NULL)
	}
	
	biopaxxml = XML::xmlInternalTreeParse(file)
	ret = list(namespaces = XML::xmlNamespaceDefinitions(XML::xmlRoot(biopaxxml), recursive = TRUE, simplify = TRUE))
	class(ret) <- c("biopax",class(ret))
	
	# deal with namespaces, we're interested in owl, rdf, biopax
	ret$ns_rdf = names(grep("rdf-syntax",ret$namespaces,ignore.case=TRUE, value=TRUE))
	ret$ns_owl = names(grep("/owl#",ret$namespaces,ignore.case=TRUE, value=TRUE))
	ret$ns_bp = names(grep("biopax-level",ret$namespaces,ignore.case=TRUE, value=TRUE))
	ret$file = file
	
	if(any(grepl("biopax-level1",ret$namespaces,ignore.case=TRUE))) {
		if(verbose) message("Found a BioPAX level 1 OWL. Unfortunatly this is not supported.\n")
		ret$biopaxlevel = 1
	}
	if(any(grepl("biopax-level2",ret$namespaces,ignore.case=TRUE))) {
		if(verbose) message("Found a BioPAX level 2 OWL. Parsing...\n")
		ret$biopaxlevel = 2
		ret$dt = internal_getBiopaxModelAsDataFrame(ret, biopaxxml, verbose=verbose)
	}
	if(any(grepl("biopax-level3",ret$namespaces,ignore.case=TRUE))) {
		if(verbose) message("Found a BioPAX level 3 OWL. Parsing...\n")
		ret$biopaxlevel = 3
		ret$dt = internal_getBiopaxModelAsDataFrame(ret, biopaxxml, verbose=verbose)
	}
	
	ret
}

#' Print a biopax object.
#'
#' @param x A \code{biopax} object to print.
#' @param ... Other arguments to be passed to \code{print}.
#' @export
#' @import data.table
#' @method print biopax
#' @examples
#'  data(biopaxexample)
#'  print(biopax) 
print.biopax <- function(x, ...) {
	
	cat("Summary of the biopax object:\n")
	print(summary(x))
	
	cat("\nInternal data:\n")
	print(x[!(names(x) %in% c("biopaxxml","df"))],...)
	
	cat(paste("Dimension of internal data.table: ",paste(dim(x$dt), collapse=","),"\n\n"))
	cat("Summary of parsed internal data.table \n")
	print(summary(x$dt))
	
	invisible(x)
}


#' This internal function parses the Biopax XML of the supplied biopax model and returns it in the data.frame format.
#' 
#' This internal function parses the Biopax XML of the supplied biopax model and returns it in the data.frame format.
#' 
#' @param biopax A biopax object 
#' @param biopaxxml Biopax XML file read in. See parseBiopax
#' @param verbose logical 
#' @return Returns the parsed biopax model in the internal data.frame format.
#' @author Frank Kramer
#' @import data.table
internal_getBiopaxModelAsDataFrame <- function (biopax, biopaxxml, verbose=TRUE) {
	
	### THIS FUNCTION ONLY RETURNS INSTANCES OF THE BIOPAX NAMESPACE! NAMESPACE MARKERS ARE STRIPPED AT THE END!
	
	## retrieve all instances below rdf:RDF 
	model = XML::getNodeSet(biopaxxml, paste("/",biopax$ns_rdf,":RDF/*",sep=""))
	
	## data.frame will be put together in the end by vectors
	if(verbose)	message("[Info Verbose] Parsing Biopax-Model as a data.table...")
	nodecount = sum(XML::xmlElementSummary(biopax$file)$nodeCounts)
	if(verbose)	message(paste("[Info Verbose] Estimating up to", nodecount,"entries. This will roughly need", round(nodecount*58*2*2/2**20), "MB of RAM."))
	if(verbose) message(paste("[Info Verbose] Where I came from this would've taken at least", round(3*nodecount/1000),"seconds!"))
	time_start = proc.time()[1]
	
	ret_colnames = c("class","id","property","property_attr","property_attr_value","property_value")
	ret = data.table(rowcount=1:nodecount,class="",id="",property="",property_attr="",property_attr_value="",property_value="", key="rowcount")
	## for each instance get data.frame entries and add them together in a df
	rowcount = 1
	for(i in 1:XML::xmlSize(model)) {
		
		## TODO Go on here: extra function. construct data.table, merge the data.tables. how to pass the dts?
		
		#class and id are name and attr
		class = XML::xmlName(model[[i]], full=T)
		id   = XML::xmlAttrs(model[[i]])[[1]]	
		#every property of an instance is represented as 1 entry in the df
		for(p in 1:XML::xmlSize(model[[i]])) {
			tryCatch({ 
						child = XML::xmlChildren(model[[i]])[[p]] 
					},
				error = function(e) { message(paste("Debug: i=",i," p=",p," class:",class, " instance:", id)) } 
			)
			
			if( (XML::xmlSize(child) > 0) && !any(class(XML::xmlChildren(child)[[1]]) %in% c("XMLInternalTextNode","XMLTextNode"))) {
				#found instancianted class here. make it a real instance, give it an id and reference it here using rdf:resource
				
				newInstanceResult = internal_XMLInstance2DF(XML::xmlChildren(child)[[1]], namespace_rdf=biopax$ns_rdf, ret, rowcount)
				if(is.null(newInstanceResult) || is.null(newInstanceResult$id)) {
					warning(paste("Something went wrong in internal_XMLInstance2DF: rowcount:", rowcount, " class:",class, " id:", id,". Debug: i=",i," p=",p))
				}
				rowcount = newInstanceResult$rowcount
				
				#generate new row to be added: reference to new instance
				tryCatch( {
					row = c(
							class, id,
							XML::xmlName(child, full=T), #property				
							paste(biopax$ns_rdf,":resource",sep=""), #property_attr
							paste("#",newInstanceResult$id,sep=""), #property_attr_value
							"" #property_value
					)
				},
				error = function(e) { 
					message(paste("Error: internal_getBioPaxModelAsDataFrame - generate new instance row:\n","p=",p," class:",class," instance:", id, " rowcount:",rowcount," newInstance:",newInstanceResult$id)) }
				)

			} else {
				attrs1 = XML::xmlAttrs(child)
				attrs2 = attributes(attrs1)
				attribute = paste(	attrs2$namespaces[[1]], attrs2$names[[1]], sep=":") #property_attr
				if(is.null(attribute)) attribute = "undefined"
				row = c(
						class, id,
						XML::xmlName(child, full=T), #property				
						attribute, #property_attr
						attrs1[[1]], #property_attr_value
						XML::xmlValue(child) #property_value
				)
			}
			
			#naive error check + add row
			if(length(row) != 6) {
				warning(paste("Something went wrong parsing: rowcount: ", rowcount, "class: ",class, "id: ", id," Parsed row: ", paste(row,collapse="|"),". Debug: i=",i," p=",p))
				
			} else {
				
				set(ret, as.integer(rowcount), 2L, row[1] )  #class
				set(ret, as.integer(rowcount), 3L, row[2] )  #id
				set(ret, as.integer(rowcount), 4L, row[3] )  #property
				set(ret, as.integer(rowcount), 5L, row[4] )  #property_attr
				set(ret, as.integer(rowcount), 6L, row[5] )  #property_attr_value
				set(ret, as.integer(rowcount), 7L, row[6] )  #property_value

				rowcount = rowcount + 1
				if(verbose)	if(rowcount%%8192 == 0) message(paste("[Info Verbose] Internal Rowcount: ",rowcount," Instance:",id,sep=" "))
			}
			# done with property
		}
		#done with instance
		#if(rowcount > 1000) break
	}

	# return result. leave out all entries that are either na or are of a different namespace than bp
	x = rowcount
	ret = copy(ret[rowcount<x & id!="",2:7, with=F])
	setkey(ret, id,class,property)
	
	# strip namespace from df$property & df$class
	set(ret, NULL, "property", stripns(ret$property))
	set(ret, NULL, "class", stripns(ret$class))
	
	if(verbose)	message(paste("[Info Verbose] Finished! Created a data.frame with ", rowcount," rows within only ",(proc.time()[1]-time_start)," seconds."))
	class(ret) <- c("biopax_df",class(ret))
	ret
}

#' This function is an internal function to count the Number of nodes and child nodes of an XMLNode.
#' 
#' This function is an internal function to count the Number of nodes and child nodes of an XMLNode.
#'  
#' @param myXMLNode XMLNode to analyze
#' @return This function returns the number of Nodes and child Nodes an XMLNode has.
#' @author Frank Kramer
#' @import data.table
internal_NrOfXMLNodes <- function(myXMLNode) {
	ret = 1
	if((XML::xmlSize(myXMLNode) > 0)  && !any(class(XML::xmlChildren(myXMLNode)[[1]]) %in% c("XMLInternalTextNode","XMLTextNode"))) {
		for(p in 1:XML::xmlSize(myXMLNode)) {
			ret = ret + internal_NrOfXMLNodes(XML::xmlChildren(myXMLNode)[[p]])
		}
	}
	ret
}

#' This function is an internal function that parses a Biopax XMLNode.
#' 
#' This function is an internal function that parses a Biopax XMLNode. Do not call it manually.
#' 
#' @param myXMLNode XMLNode
#' @param namespace_rdf String specifying the namespace to use for rdf:resource and rdf:datatype
#' @param ret data.table object contaning the already parsed data to attach this instance to
#' @param rowcount Numeric specifying the row at which further parsed data is inserted into the data.table
#' @return Returns a list contianing the new rowcount and the instance id of the added instance
#' @author Frank Kramer
#' @import data.table
internal_XMLInstance2DF <- function(myXMLNode, namespace_rdf, ret, rowcount) {

	class = XML::xmlName(myXMLNode, full=T)
	id   = XML::xmlAttrs(myXMLNode)[[1]]	
	#every property of an instance is represented as 1 entry in the df
	
	for(x in 1:XML::xmlSize(myXMLNode)) {
		tryCatch({ 
					child = XML::xmlChildren(myXMLNode)[[x]] 
				},
				error = function(e) { message(paste("Debug: internal_XMLInstance2DF xmlChildren x=",x," class:",class, " instance:", id)) } 
		)
		if( (XML::xmlSize(child) > 0) && !any(class(XML::xmlChildren(child)[[1]]) %in% c("XMLInternalTextNode","XMLTextNode"))) {
			
			#found instancianted class here. make it a real instance, give it an id and reference it here using rdf:resource
			newInstanceResult = internal_XMLInstance2DF(XML::xmlChildren(child)[[1]], namespace_rdf=namespace_rdf, ret, rowcount)
			if(is.null(newInstanceResult) || is.null(newInstanceResult$id)) {
				warning(paste("Something went wrong in internal_XMLInstance2DF: rowcount:", rowcount, " class:",class, " id:", id,". Debug: i=",i," p=",p))
			}
			rowcount = newInstanceResult$rowcount
			
			
			#generate new row to be added: reference to new instance
			tryCatch( {
						row = c(
								class, id,
								XML::xmlName(child, full=T), #property				
								paste(namespace_rdf,":resource",sep=""), #property_attr
								paste("#",newInstanceResult$id,sep=""), #property_attr_value
								"" #property_value
						)
					},
					error = function(e) { 
						message(paste("Error: internal_XMLInstance2DF - generate new instance row:\n","x=",x," class:",class," instance:", id, " rowcount:",rowcount," newInstance:",newInstanceResult$id)) }
			)
			
		}
		else {
			attrs1 = XML::xmlAttrs(child)
			attrs2 = attributes(attrs1)
			attribute = paste(	attrs2$namespaces[[1]], attrs2$names[[1]], sep=":") #property_attr
			if(is.null(attribute)) attribute = "undefined"
			row = c(
					class, id,
					XML::xmlName(child, full=T), #property				
					attribute, #property_attr
					attrs1[[1]], #property_attr_value
					XML::xmlValue(child) #property_value
			)
		}
		
		#naive error check + add row
		if(length(row) != 6) {
			warning(paste("Something went wrong parsing: rowcount:", rowcount, "class:",class, "id:", id," Parsed row:", paste(row,collapse="|"),". Debug: x=",x))
		} else {
			
			set(ret, as.integer(rowcount), 2L, row[1] )  #class
			set(ret, as.integer(rowcount), 3L, row[2] )  #id
			set(ret, as.integer(rowcount), 4L, row[3] )  #property
			set(ret, as.integer(rowcount), 5L, row[4] )  #property_attr
			set(ret, as.integer(rowcount), 6L, row[5] )  #property_attr_value
			set(ret, as.integer(rowcount), 7L, row[6] )  #property_value
			
			rowcount = rowcount + 1
		}
		# done with child
	}
	list(rowcount=rowcount, id=id)
}
################### PARSING END ###################
