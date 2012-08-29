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
#' 		\item{biopaxxml}{NULL}
#' 		\item{summary}{NULL}
#' 		\item{df}{The data.frame representing the biopax in R}
#' 		\item{ns_rdf}{RDF Namespace}
#'		\item{ns_owl}{OWL Namespace}
#'		\item{ns_bp}{Biopax Namespace}
#' 		\item{file}{NULL}
#'	}
#'  
#' @param level integer. Specifies the BioPAX level. 
#' @returnType A biopax model
#' @return A biopax model
#' @author Frank Kramer
#' @export
createBiopax <- function(level = 2)  {
	ret = list(biopaxxml = NULL)
	ret$summary = NULL
	ret$file = NULL
	class(ret) <- c("biopax",class(ret))
	
	if(level==1) cat("BioPAX level 1 OWL is unfortunatly not supported.\n")
	if(level==2) {
		ret$biopaxlevel = 2
		df_colnames = c("instancetype","instanceid","property","property_attr","property_attr_value","property_value")
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
	if(level==3) cat("BioPAX level 3 OWL is not yet supported but work is underway.\n")
	
	ret
}


#' This function reads in a Biopax .owl file
#' 
#' This function reads in a Biopax .owl file and generates the internal data.frame format used in this package.
#' This function can take a while with really big Biopax files like NCIs Pathway Interaction Database or Reactome.
#' In almost every case this is your starting point.
#' Returns a biopax model, which is a list with named elements: 
#' \describe{
#' 		\item{biopaxxml}{The XML which was read in from file.}
#' 		\item{summary}{The generated summary statistic.}
#' 		\item{df}{The data.frame representing the biopax in R}
#' 		\item{ns_rdf}{RDF Namespace}
#'		\item{ns_owl}{OWL Namespace}
#'		\item{ns_bp}{Biopax Namespace}
#' 		\item{file}{File name}
#'	}
#'  
#' @param file string. File name 
#' @param verbose logical. Output messages about how parsing is going and so on.
#' @param generateSummary logical. Generates and attaches the summary for the Biopax file to biopax$summary
#' @param generateDF logical. If this is set to FALSE no data.frame is generated, only the file is read in to biopax$biopaxxml
#' @returnType A biopax model
#' @return A biopax model
#' @author Frank Kramer
#' @export
readBiopax <- function(file, verbose=TRUE, generateSummary=TRUE, generateDF=TRUE)  {
	ret = list(biopaxxml = XML::xmlInternalTreeParse(file))
	class(ret) <- c("biopax",class(ret))
	ret$namespaces = XML::xmlNamespaceDefinitions(XML::xmlRoot(ret$biopaxxml), recursive = TRUE, simplify = TRUE)
	
	# deal with namespaces, we're interested in owl, rdf, biopax
	ret$ns_rdf = names(grep("rdf-syntax",ret$namespaces,ignore.case=TRUE, value=TRUE))
	ret$ns_owl = names(grep("/owl#",ret$namespaces,ignore.case=TRUE, value=TRUE))
	ret$ns_bp = names(grep("biopax-level",ret$namespaces,ignore.case=TRUE, value=TRUE))
	ret$file = file
	
	if(any(grepl("biopax-level1",ret$namespaces,ignore.case=TRUE))) {
		if(verbose) cat("Found a BioPAX level 1 OWL. Unfortunatly this is not supported.\n")
		ret$biopaxlevel = 1
		if(generateSummary)	ret$summary = summary_biopax(ret, verbose=verbose)
	}
	if(any(grepl("biopax-level2",ret$namespaces,ignore.case=TRUE))) {
		if(verbose) cat("Found a BioPAX level 2 OWL. Parsing...\n")
		ret$biopaxlevel = 2
		if(generateSummary)	ret$summary = summary_biopax(ret, verbose=verbose)
		if(generateDF) ret$df = internal_getBiopaxModelAsDataFrame(ret, verbose=verbose)
	}
	if(any(grepl("biopax-level3",ret$namespaces,ignore.case=TRUE))) {
		if(verbose) cat("Found a BioPAX level 3 OWL. Unfortunatly this is not supported yet, but support will be added in the near future!\n")
		ret$biopaxlevel = 3
		if(generateSummary)	ret$summary = summary_biopax(ret, verbose=verbose)
	}
	
	ret
}

#' This function generates a summary statistic for the biopax model.
#' 
#' This function generates a summary statistic for the biopax model.
#' This function can take a while with really big Biopax files like NCIs Pathway Interaction Database or Reactome.
#' This function is called internally by readBiopax if generateSummary == TRUE, so just check biopax$summary if it has already been generated. 
#' 
#' @param object A biopax model 
#' @param verbose logical
#' @returnType biopax2_summary
#' @return Returns the summary for the supplied biopax model.
#' @author Frank Kramer
#' @export
summary_biopax <- function(object, verbose=TRUE) {
	model = XML::getNodeSet(object$biopaxxml,"/rdf:RDF/*")
	modellength = length(model)
	if(verbose)	cat("Generating summary for model with ", modellength,"entries.\n")
	nodenames = list()
	for(i in 1:modellength) {
		nodenames[i] = XML::xmlName(model[[i]], full=T)
	}	
	
	nodenames = factor(unlist(nodenames))
	nodenames_factors = levels(nodenames)
	nodenames_summary = data.frame(instancetype = nodenames_factors, count=NA)
	
	for(i in 1:length(nodenames_factors)) {
		nodenames_summary[i,2] = sum(nodenames==nodenames_factors[i])
	}
	if(verbose)	cat("Done!\n")
	ret = list(size=modellength,instances=nodenames_summary)
	class(ret) <- c("biopax_summary",class(ret))
	ret
}

#' This function returns a summary statistic for the biopax model.
#' 
#' This function returns a summary statistic for the biopax model. 
#' It checks if the summary has already been generated with summary.biopax and generates it if necessary.
#' Have a look at the biopax$summary entry for more statistics.
#' This function can take a while with really big Biopax files like NCIs Pathway Interaction Database or Reactome.
#' 
#' @param biopax A biopax model 
#' @returnType summary 
#' @return This function returns a summary statistic for the biopax model.
#' @author Frank Kramer
#' @export
listBiopaxNodes <- function(biopax) {
	if(is.null(biopax$summary)) {
		biopax$summary = summary_biopax(biopax,verbose=TRUE)
	}
	biopax$summary$instances
}

#' This internal function parses the Biopax XML of the supplied biopax model and returns it in the data.frame format.
#' 
#' This internal function parses the Biopax XML of the supplied biopax model and returns it in the data.frame format.
#' 
#' @param biopax A biopax model 
#' @param verbose logical 
#' @returnType data.frame
#' @return Returns the parsed biopax model in the internal data.frame format.
#' @author Frank Kramer
#  #' @export
internal_getBiopaxModelAsDataFrame <- function (biopax, verbose=TRUE) {
	
	### THIS FUNCTION ONLY RETURNS INSTANCES OF THE BIOPAX NAMESPACE! NAMESPACE MARKERS ARE STRIPPED AT THE END!
	
	## retrieve all instances below rdf:RDF 
	model = XML::getNodeSet(biopax$biopaxxml, paste("/",biopax$ns_rdf,":RDF/*",sep=""))
	
	## data.frame will be put together in the end by vectors
	if(verbose)	cat("Parsing Biopax-Model as a data.frame...\n")
	nodecount = sum(XML::xmlElementSummary(biopax$file)$nodeCounts)
	if(verbose)	cat("Estimating up to", nodecount,"entries. \nThis will roughly need", round(nodecount*58*2*2/2**20), "MB of RAM.\n")
	if(verbose)	cat("Where I came from this would've taken at least", round(nodecount/1000),"seconds!\n")
	time_start = proc.time()[1]
	ret = matrix(data=NA,nrow = nodecount, ncol = 6)
	
	ret_colnames = c("instancetype","instanceid","property","property_attr","property_attr_value","property_value")
	colnames(ret) = ret_colnames
	## for each instance get data.frame entries and add them together in a df
	rowcount = 0
	for(i in 1:XML::xmlSize(model)) {
		
		#instancetype and id are name and attr
		instancetype = XML::xmlName(model[[i]], full=T)
		instanceid   = XML::xmlAttrs(model[[i]])[[1]]	
		#every property of an instance is represented as 1 entry in the df
		for(p in 1:XML::xmlSize(model[[i]])) {
			tryCatch({ 
						child = XML::xmlChildren(model[[i]])[[p]] 
					},
				error = function(e) { print(paste("Debug: i=",i," p=",p," instance:",instancetype, instanceid)) } 
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
							print(paste("Error: internal_getBioPaxModelAsDataFrame - NewInstance:","x=",x,"p=",p," instance:",instancetype, instanceid));
							print("ret:"); print(ret[1:(max(1,rowcount)),]); 
							print("newInstance:"); print(newInstance); print(rowcount) }
					)
				}

				#generate new row to be added: reference to new instance
				tryCatch( {
					row = c(
							instancetype, instanceid,
							XML::xmlName(child, full=T), #property				
							paste(biopax$ns_rdf,":resource",sep=""), #property_attr
							paste("#",newInstance[1,"instanceid"],sep=""), #property_attr_value
							"" #property_value
					)
				},
				error = function(e) { 
					print(paste("Error: internal_getBioPaxModelAsDataFrame - NewInstance:","p=",p," instance:",instancetype, instanceid));
					print("ret:"); print(ret[1:(max(1,rowcount)),]); 
					print("newInstance:"); print(newInstance); print(rowcount) }
				)

			} else {
				row = c(
						instancetype, instanceid,
						XML::xmlName(child, full=T), #property				
						paste(	attributes(XML::xmlAttrs(child))$namespaces[[1]],
								attributes(XML::xmlAttrs(child))$names[[1]],
								sep=":"), #property_attr
						XML::xmlAttrs(child)[[1]], #property_attr_value
						XML::xmlValue(child) #property_value
				)
			}
			
			#naive error check + add row
			if(length(row) != 6) {
				warning(paste("Something went wrong parsing:",instancetype, instanceid," Parsed ",row,". Debug: i=",i," p=",p))
				
			} else {
				rowcount = rowcount + 1
				ret[rowcount,] = row
				if(verbose)	if(rowcount%%8192 == 0) cat(paste("Rowcount: ",rowcount,"InstanceNr:",i,"Instance:",instanceid,"\n",sep=" "))
				#if(verbose)	if(rowcount%%100 == 0) cat(paste("Rowcount: ",rowcount,"InstanceNr:",i,"Instance:",instanceid,"\n",sep=" "))
			}
			# done with property
		}
		#done with instance
		#if(rowcount > 1000) break
	}
	# return result. leave out all entries that are either na or are of a different namespace than bp
	ret2 = data.frame( ret[ !(is.na(ret[,1]) | !(isOfNamespace(ret[,1],biopax$ns_bp)) ) ,] )
	# strip namespace from df$property & df$instancetype
	levels(ret2$property) =  stripns(levels(ret2$property))
	levels(ret2$instancetype) =  stripns(levels(ret2$instancetype))

	# take care if matrix has only one row
	if(dim(ret2)[[1]]==6 & dim(ret2)[[2]]==1) {
		ret2=t(ret2)
		rownames(ret2)=1
	}

	colnames(ret2) = ret_colnames
	rm(ret)
	if(verbose)	cat(paste("Finished! Created a data.frame with", rowcount,"rows within only",(proc.time()[1]-time_start),"seconds.\n"))
	class(ret2) <- c("biopax2_df",class(ret2))
	ret2
}

#' This function in an internal function to count the Number of nodes and child nodes of an XMLNode.
#' 
#' This function in an internal function to count the Number of nodes and child nodes of an XMLNode.
#'  
#' @param myXMLNode XMLNode to analyze
#' @returnType integer
#' @return This function returns the number of Nodes and child Nodes an XMLNode has.
#' @author Frank Kramer
#  #' @export
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
#' @returnType matrix
#' @return Returns the matrix generated by parsing the XMLNode
#' @author Frank Kramer
#  #' @export
internal_XMLInstance2DF <- function(myXMLNode, namespace_rdf) {
	ret_colnames = c("instancetype","instanceid","property","property_attr","property_attr_value","property_value")
	ret = matrix(data=NA,nrow = 1000, ncol = 6)
	colnames(ret) = ret_colnames
	instancetype = XML::xmlName(myXMLNode, full=T)
	instanceid   = XML::xmlAttrs(myXMLNode)[[1]]	
	#every property of an instance is represented as 1 entry in the df
	rowcount = 0
	for(x in 1:XML::xmlSize(myXMLNode)) {
		tryCatch({ 
					child = XML::xmlChildren(myXMLNode)[[x]] 
				},
				error = function(e) { print(paste("Debug: internal_XMLInstance2DF xmlChildren x=",x," instance:",instancetype, instanceid)) } 
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
							print(paste("Error: internal_XMLInstance2DF - NewInstance:","x=",x,"y=",y," instance:",instancetype, instanceid));
							print("ret:"); print(ret[1:(max(1,rowcount)),]); 
							print("newInstance:"); print(newInstance); print(rowcount) }
				)
			}
			#generate new row to be added: reference to new instance
			tryCatch( {
						row = c(
								instancetype, instanceid,
								XML::xmlName(child, full=T), #property				
								paste(namespace_rdf,":resource",sep=""), #property_attr
								paste("#",newInstance[1,"instanceid"],sep=""), #property_attr_value
								"" #property_value
						)
					},
					error = function(e) { 
						print(paste("Error: internal_XMLInstance2DF - NewInstance:","x=",x," instance:",instancetype, instanceid, XML::xmlName(child, full=T)));
						print("ret:"); print(ret[1:(max(1,rowcount)),]); 
						print("newInstance:"); print(newInstance); print(rowcount) }
			)
			
		}
		else {
			row = c(
					instancetype, instanceid,
					XML::xmlName(child, full=T), #property				
					paste(	attributes(XML::xmlAttrs(child))$namespaces[[1]],
							attributes(XML::xmlAttrs(child))$names[[1]],
							sep=":"), #property_attr
					XML::xmlAttrs(child)[[1]], #property_attr_value
					XML::xmlValue(child) #property_value
			)
		}
		
		#naive error check + add row
		if(length(row) != 6) {
			warning(paste("Something went wrong parsing:",instancetype, instanceid,". Parsed ",row,". It was an instance created inside of a property. Debug: x=",x," y=",y))
		} else {
			rowcount = rowcount + 1
			ret[rowcount,] = row
		}
		# done with child
	}
	matrix( ret[!is.na(ret[,1]),], ncol=6, dimnames=list(list(),ret_colnames))
}
################### PARSING END ###################
