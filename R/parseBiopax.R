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

##TODO: FIX PARSING. SEE PATHWAY COMMONS

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
#' @examples
#'  biopax = createBiopax(level=2)
#' 
createBiopax <- function(level = 2)  {
	ret = list(file = NULL)
	class(ret) <- c("biopax",class(ret))
	
	if(level==1) message("BioPAX level 1 OWL is unfortunatly not supported.")
	if(level==2) {
		ret$biopaxlevel = 2
		df_colnames = c("class","id","property","property_attr","property_attr_value","property_value")
		ret$df = data.frame(rbind(df_colnames,0))
		colnames(ret$df) = df_colnames
		ret$df = ret$df[0,]

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
		df_colnames = c("class","id","property","property_attr","property_attr_value","property_value")
		ret$df = data.frame(rbind(df_colnames,0))
		colnames(ret$df) = df_colnames
		ret$df = ret$df[0,]
		
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
		ret$df = internal_getBiopaxModelAsDataFrame(ret, biopaxxml, verbose=verbose)
	}
	if(any(grepl("biopax-level3",ret$namespaces,ignore.case=TRUE))) {
		if(verbose) message("Found a BioPAX level 3 OWL. Parsing...\n")
		ret$biopaxlevel = 3
		ret$df = internal_getBiopaxModelAsDataFrame(ret, biopaxxml, verbose=verbose)
	}
	
	ret
}

#' Print a biopax object.
#'
#' @param x A \code{biopax} object to print.
#' @param ... Other arguments to be passed to \code{print}.
#' @export
#' @method print biopax
#' @examples
#'  data(biopax2example)
#'  print(biopax) 
print.biopax <- function(x, ...) {
	
	cat("Summary of the biopax object:\n")
	print(summary(x))
	
	cat("\nInternal data:\n")
	print(x[!(names(x) %in% c("biopaxxml","df"))],...)
	
	cat(paste("Dimension of internal data.frame: ",paste(dim(x$df), collapse=","),"\n\n"))
	cat("Summary of parsed internal data.frame \n")
	print(summary(x$df))
	
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
internal_getBiopaxModelAsDataFrame <- function (biopax, biopaxxml, verbose=TRUE) {
	
	### THIS FUNCTION ONLY RETURNS INSTANCES OF THE BIOPAX NAMESPACE! NAMESPACE MARKERS ARE STRIPPED AT THE END!
	
	## retrieve all instances below rdf:RDF 
	model = XML::getNodeSet(biopaxxml, paste("/",biopax$ns_rdf,":RDF/*",sep=""))
	
	## data.frame will be put together in the end by vectors
	if(verbose)	message("[Info Verbose] Parsing Biopax-Model as a data.frame...")
	nodecount = sum(XML::xmlElementSummary(biopax$file)$nodeCounts)
	if(verbose)	message(paste("[Info Verbose] Estimating up to", nodecount,"entries. This will roughly need", round(nodecount*58*2*2/2**20), "MB of RAM."))
	if(verbose) message(paste("[Info Verbose] Where I came from this would've taken at least", round(3*nodecount/1000),"seconds!"))
	time_start = proc.time()[1]
	ret = matrix(data=NA,nrow = nodecount, ncol = 6)
	
	ret_colnames = c("class","id","property","property_attr","property_attr_value","property_value")
	colnames(ret) = ret_colnames
	## for each instance get data.frame entries and add them together in a df
	rowcount = 0
	for(i in 1:XML::xmlSize(model)) {
		
		#class and id are name and attr
		class = XML::xmlName(model[[i]], full=T)
		id   = XML::xmlAttrs(model[[i]])[[1]]	
		#every property of an instance is represented as 1 entry in the df
		for(p in 1:XML::xmlSize(model[[i]])) {
			tryCatch({ 
						child = XML::xmlChildren(model[[i]])[[p]] 
					},
				error = function(e) { message(paste("Debug: i=",i," p=",p," instance:",class, id)) } 
			)
			
			if( (XML::xmlSize(child) > 0) && !any(class(XML::xmlChildren(child)[[1]]) %in% c("XMLInternalTextNode","XMLTextNode"))) {
				#found instancianted class here. make it a real instance, give it an id and reference it here using rdf:resource
				# can check with i = 7290, p=3 on PID
				newInstance = internal_XMLInstance2DF(XML::xmlChildren(child)[[1]], namespace_rdf=biopax$ns_rdf)
				for(x in 1:dim(newInstance)[[1]]) {
					rowcount = rowcount + 1
					tryCatch( {
						ret[rowcount,] = matrix(newInstance,ncol=6)[x,]
						},
						error = function(e) { 
							message(paste("Error: internal_getBioPaxModelAsDataFrame - NewInstance:","x=",x,"p=",p," instance:",class, id));
							message(paste("ret:",ret[1:(max(1,rowcount)),],"newInstance:",newInstance, rowcount))}
					)
				}

				#generate new row to be added: reference to new instance
				tryCatch( {
					row = c(
							class, id,
							XML::xmlName(child, full=T), #property				
							paste(biopax$ns_rdf,":resource",sep=""), #property_attr
							paste("#",newInstance[1,"id"],sep=""), #property_attr_value
							"" #property_value
					)
				},
				error = function(e) { 
					message(paste("Error: internal_getBioPaxModelAsDataFrame - NewInstance:","p=",p," instance:",class, id));
					message(paste("ret:",ret[1:(max(1,rowcount)),],"newInstance:",newInstance,rowcount)) }
				)

			} else {
				attribute = paste(	attributes(XML::xmlAttrs(child))$namespaces[[1]], attributes(XML::xmlAttrs(child))$names[[1]], sep=":") #property_attr
				if(is.null(attribute)) attribute = "undefined"
				row = c(
						class, id,
						XML::xmlName(child, full=T), #property				
						attribute, #property_attr
						XML::xmlAttrs(child)[[1]], #property_attr_value
						XML::xmlValue(child) #property_value
				)
			}
			
			#naive error check + add row
			if(length(row) != 6) {
				warning(paste("Something went wrong parsing:",class, id," Parsed ",row,". Debug: i=",i," p=",p))
				
			} else {
				rowcount = rowcount + 1
				ret[rowcount,] = row
				if(verbose)	if(rowcount%%8192 == 0) message(paste("[Info Verbose] Internal Rowcount: ",rowcount,"InstanceNr:",i,"Instance:",id,sep=" "))
			}
			# done with property
		}
		#done with instance
		#if(rowcount > 1000) break
	}
	# return result. leave out all entries that are either na or are of a different namespace than bp
	ret2 = data.frame( ret[ !(is.na(ret[,1]) | !(isOfNamespace(ret[,1],biopax$ns_bp)) ) ,] )
	# strip namespace from df$property & df$class
	levels(ret2$property) =  stripns(levels(ret2$property))
	levels(ret2$class) =  stripns(levels(ret2$class))

	# take care if matrix has only one row
	if(dim(ret2)[[1]]==6 & dim(ret2)[[2]]==1) {
		ret2=t(ret2)
		rownames(ret2)=1
	}

	colnames(ret2) = ret_colnames
	rm(ret)
	if(verbose)	message(paste("[Info Verbose] Finished! Created a data.frame with", rowcount,"rows within only",(proc.time()[1]-time_start),"seconds."))
	class(ret2) <- c("biopax_df",class(ret2))
	ret2
}

#' This function in an internal function to count the Number of nodes and child nodes of an XMLNode.
#' 
#' This function in an internal function to count the Number of nodes and child nodes of an XMLNode.
#'  
#' @param myXMLNode XMLNode to analyze
#' @return This function returns the number of Nodes and child Nodes an XMLNode has.
#' @author Frank Kramer
internal_NrOfXMLNodes <- function(myXMLNode) {
	ret = 1
	if((XML::xmlSize(myXMLNode) > 0)  && !any(class(XML::xmlChildren(myXMLNode)[[1]]) %in% c("XMLInternalTextNode","XMLTextNode"))) {
		for(p in 1:XML::xmlSize(myXMLNode)) {
			ret = ret + internal_NrOfXMLNodes(XML::xmlChildren(myXMLNode)[[p]])
		}
	}
	ret
}

#' This function in an internal function that parses a Biopax Level 2 XMLNode.
#' 
#' This function in an internal function that parses a Biopax Level 2 XMLNode.
#' 
#' @param myXMLNode XMLNode
#' @param namespace_rdf String specifying the namespace to use for rdf:resource and rdf:datatype
#' @return Returns the matrix generated by parsing the XMLNode
#' @author Frank Kramer
internal_XMLInstance2DF <- function(myXMLNode, namespace_rdf) {
	ret_colnames = c("class","id","property","property_attr","property_attr_value","property_value")
	ret = matrix(data=NA,nrow = 10000, ncol = 6)
	colnames(ret) = ret_colnames
	class = XML::xmlName(myXMLNode, full=T)
	id   = XML::xmlAttrs(myXMLNode)[[1]]	
	#every property of an instance is represented as 1 entry in the df
	rowcount = 0
	for(x in 1:XML::xmlSize(myXMLNode)) {
		tryCatch({ 
					child = XML::xmlChildren(myXMLNode)[[x]] 
				},
				error = function(e) { message(paste("Debug: internal_XMLInstance2DF xmlChildren x=",x," instance:",class, id)) } 
		)
		if( (XML::xmlSize(child) > 0) && !any(class(XML::xmlChildren(child)[[1]]) %in% c("XMLInternalTextNode","XMLTextNode"))) {
			#found instancianted class here. make it a real instance, give it an id and referece it here using rdf:resource
			newInstance = internal_XMLInstance2DF(XML::xmlChildren(child)[[1]])
			for(y in 1:dim(newInstance)[[1]]) {
				rowcount = rowcount + 1
				tryCatch( {
							ret[rowcount,] = matrix(newInstance,ncol=6)[y,]
						},
						error = function(e) { 
							message(paste("Error: internal_XMLInstance2DF - NewInstance:","x=",x,"y=",y," instance:",class, id));
							message(paste("ret:",ret[1:(max(1,rowcount)),],"newInstance:",newInstance,rowcount)) }
				)
			}
			#generate new row to be added: reference to new instance
			tryCatch( {
						row = c(
								class, id,
								XML::xmlName(child, full=T), #property				
								paste(namespace_rdf,":resource",sep=""), #property_attr
								paste("#",newInstance[1,"id"],sep=""), #property_attr_value
								"" #property_value
						)
					},
					error = function(e) { 
						message(paste("Error: internal_XMLInstance2DF - NewInstance:","x=",x," instance:",class, id, XML::xmlName(child, full=T)));
						message(paste(("ret:"),ret[1:(max(1,rowcount)),],"newInstance:",newInstance,rowcount)) }
			)
			
		}
		else {
			attribute = paste(	attributes(XML::xmlAttrs(child))$namespaces[[1]], attributes(XML::xmlAttrs(child))$names[[1]], sep=":") #property_attr
			if(is.null(attribute)) attribute = "undefined"
			row = c(
					class, id,
					XML::xmlName(child, full=T), #property				
					attribute, #property_attr
					XML::xmlAttrs(child)[[1]], #property_attr_value
					XML::xmlValue(child) #property_value
			)
		}
		
		#naive error check + add row
		if(length(row) != 6) {
			warning(paste("Something went wrong parsing:",class, id,". Parsed ",row,". It was an instance created inside of a property. Debug: x=",x," y=",y))
		} else {
			rowcount = rowcount + 1
			ret[rowcount,] = row
		}
		# done with child
	}
	matrix( ret[!is.na(ret[,1]),], ncol=6, dimnames=list(list(),ret_colnames))
}
################### PARSING END ###################
