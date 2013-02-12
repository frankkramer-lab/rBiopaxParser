###############################################################################
#
# visualizeBiopax.R: 	This file contains the all functions related to visualizing a parsed Biopax model.
#						Generation of graphs, layout, adjacency matrices, modification of graphs.
# author: Frank Kramer <dev@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################


#' This function generates an adjacency matrix from the activations/inhibitions of a pathway in a biopax model.
#'  
#' This function internally first calls pathway2RegulatoryGraph, then converts the regulatory graph to an adjacency matrix.
#' See pathway2RegulatoryGraph for more details.
#'  
#' @param biopax A biopax model
#' @param pwid string
#' @param expandSubpathways logical. If TRUE subpathways are expanded into this graph, otherwise only this very pathway is used.
#' @param splitComplexMolecules logical. If TRUE every complex is split up into its components. This leads to splitting a single node with name of the complex into several nodes with names of the components, these components all have identical edges.   
#' @param useIDasNodenames logical. If TRUE nodes of the graph are named by their molecule IDs instead of using the NAME property. This can help with badly annotated/formatted databases.
#' @param verbose logical 
#' @return Returns the adjacency matrix representing the regulatory graph of the supplied pathway.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  pwid1 = "pid_p_100002_wntpathway"
#'  pwid2 = "pid_p_100146_hespathway"
#'  pathway2AdjacancyMatrix(biopax, pwid1)

pathway2AdjacancyMatrix <- function(biopax, pwid, expandSubpathways=TRUE, splitComplexMolecules=TRUE, useIDasNodenames=FALSE, verbose=TRUE) {
	
	if(!require(graph)) {
		message(paste("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!","\n"))
		return(NULL)
	}
	
	mygraph = pathway2RegulatoryGraph(biopax, pwid, expandSubpathways=expandSubpathways, splitComplexMolecules=splitComplexMolecules, useIDasNodenames=useIDasNodenames, verbose=verbose)
	ttt = as(mygraph,"graphAM")
	as(ttt,"matrix")
	
}

#test: pwid = "pid_p_200067_wnt_signaling_pathway"
#      i    = "pid_i_203122_inhibitor_1"
#      pwid = "pid_p_200020_wnt_noncanonical_pathway"
#      i    = "pid_i_203122_inhibitor_1"
#     

#' This function generates the regulatory graph from the activations/inhibitions of a pathway in a biopax model.
#'  
#' This functions builds a graph from the pathway components of the supplied pathway. 
#' Only instances of class 'control' are considered, this leads a functinal graph with all edges either representing activations or inhibitions. No transports, no translocation, etc.
#' If desired complexes can be split up into several nodes, this can sometimes lead to a more complex and cluttered graph.
#' There can not be multiple edges between 2 nodes. Whenever duplicated edges are generated (especially by splitting up complexes) a warning is thrown.
#'  
#' @param biopax A biopax model
#' @param pwid string
#' @param expandSubpathways logical. If TRUE subpathways are expanded into this graph, otherwise only this very pathway is used.
#' @param splitComplexMolecules logical. If TRUE every complex is split up into its components. This leads to splitting a single node with name of the complex into several nodes with names of the components, these components all have identical edges.   
#' @param useIDasNodenames logical. If TRUE nodes of the graph are named by their molecule IDs instead of using the NAME property. This can help with badly annotated/formatted databases.
#' @param verbose logical 
#' @return Returns the representing the regulatory graph of the supplied pathway in a node-edge-list graph.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  pwid1 = "pid_p_100002_wntpathway"
#'  pwid2 = "pid_p_100146_hespathway"
#'  mygraph = pathway2RegulatoryGraph(biopax, pwid1)
#'  plotRegulatoryGraph(mygraph)
pathway2RegulatoryGraph  <- function(biopax, pwid, expandSubpathways=TRUE, splitComplexMolecules=TRUE, useIDasNodenames=FALSE, verbose=TRUE) {
	
	if(!require(graph)) {
		message(paste("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!","\n"))
		return(NULL)
	}
	
	pwid = unique(striphash(pwid))
	#first make an node edge list, because its easier
	mygraph = new(getClassDef("graphNEL", package="graph"))
	graph::edgemode(mygraph) = "directed"
	
	#get pathway component list
	pw_component_list = selectInstances(biopax,id=listPathwayComponents(biopax,pwid)$id)
	#consider only controls //TODO: Consider subpathways!
	for(i in unique(pw_component_list$id)) {
		instance = pw_component_list[pw_component_list$id==i,]
		if(any(isOfClass(instance,"control", considerInheritance = TRUE))) {
			#if type neither inhibition nor activation i dont know what to do anyways. 
			#Cancer Cell Map doesnt set control-type of catalysises to ACTIVATION although it is required by biopax level 2standards. quick fix here:
			if(biopax$biopaxlevel == 2) {
				type = as.character(instance[tolower(instance$property)==tolower("CONTROL-TYPE"),"property_value"])
			}
			if(biopax$biopaxlevel == 3) {
				type = as.character(instance[tolower(instance$property)==tolower("controlType"),"property_value"])
			}
			if(length(type)==0) {
				if(any(isOfClass(instance,"catalysis"))) {
					type="ACTIVATION"
				} else {
					next
				}
			}
			if ( !grepl("activation",type,ignore.case=TRUE) & !grepl("inhibition",type,ignore.case=TRUE)) {
				next
			}
			
			#controllers are physicalentityparticipants, get those
			controller_ids = as.character(unique(instance[tolower(instance$property)==tolower("CONTROLLER"),"property_attr_value"]))
			controllers = NA
			for(i2 in controller_ids) {
				c_instance = selectInstances(biopax, id=i2)
				#each physicalentityparticipant has exactly 1 physicalentity, get that and get the name!
				if(biopax$biopaxlevel == 2) {
					c_instance = selectInstances(biopax, id=c_instance[tolower(c_instance$property)==tolower("PHYSICAL-ENTITY"),"property_attr_value"])
				}
				if(splitComplexMolecules & any(isOfClass(c_instance,"complex"))) {
					if(useIDasNodenames) {
						controllers = unique(c(controllers,  as.character( splitComplex(biopax,i2)$id )))
					} else {
						controllers = unique(c(controllers,  as.character( splitComplex(biopax,i2)$name )))						
					}
				} else {
					if(useIDasNodenames) {
						controllers = unique(c(controllers, as.character(c_instance$id[1])))
					} else {
						controllers = unique(c(controllers, getInstanceProperty(biopax,c_instance$id[1])))						
					}
				}	
			}
			
			#controlleds are interactions or pathways. ignoring pathways for now
			controlled_ids = as.character(unique(instance[tolower(instance$property)==tolower("CONTROLLED"),"property_attr_value"]))
			controlleds = NA
			for(i2 in controlled_ids) {
				c_instance = selectInstances(biopax, id=i2)
				#depending on type of c_instance we must differentiate here. 
				#pathway=ignore,interaction=ignore(this is only a nice field with a name in PID anyways),
				#complexAssembly=we deal with a complex. split up or dont and return names
				#biochemicalReaction=...
				#any conversion: add up left & rights and get the names
				if(any(isOfClass(c_instance,c("conversion"),considerInheritance = TRUE))) {
					leftrights = as.character(c_instance[tolower(c_instance$property)==tolower("LEFT") | tolower(c_instance$property)==tolower("RIGHT") ,"property_attr_value"])
					for(i3 in leftrights) {
						#every left/right is an physicalentityparticipants, get that as above
						leftrights_instance = selectInstances(biopax, id=i3)
						if(biopax$biopaxlevel == 2) {
							leftrights_instance = selectInstances(biopax, id=leftrights_instance[tolower(leftrights_instance$property)==tolower("PHYSICAL-ENTITY"),"property_attr_value"])
						}
						#split complexes?
						if(splitComplexMolecules & any(isOfClass(leftrights_instance,"complex"))) {
							if(useIDasNodenames) {
								controlleds = unique(c(controlleds,  as.character( splitComplex(biopax,i3)$id )))
							} else {
								controlleds = unique(c(controlleds,  as.character( splitComplex(biopax,i3)$name )))						
							}
						} else {
							if(useIDasNodenames) {
								controlleds = unique(c(controlleds, as.character(leftrights_instance$id[1])))
							} else {
								controlleds = unique(c(controlleds, getInstanceProperty(biopax,leftrights_instance$id[1])))						
							}
						}
					}
				}
			}
			
			controllers = controllers[!is.na(controllers) & !is.null(controllers) & nchar(controllers) > 0 ]
			controlleds = controlleds[!is.na(controlleds) & !is.null(controlleds) & nchar(controlleds) > 0 ]
			#only add viable directed edges:
			if(length(controllers)==0 | length(controlleds)==0) {
				next
			}
			
			#verbose
			if(verbose) {
				message(paste("Adding to graph: ",unique(as.character(instance[,"class"])), "-", type,
								"Controllers: ", paste(controllers, collapse=" "), "Controlleds: ", paste(controlleds, collapse=" ")))
			}
			#now we have the controllers and controlleds lists, add edge
			# leave out duplicates
			newnodes = unique(c(controllers,controlleds))
			# only add new nodes that are not null, na or empty strings
			newnodes = newnodes[!(newnodes %in% graph::nodes(mygraph))]
			if(length(newnodes)>0) mygraph = graph::addNode(newnodes,mygraph)
			#set weight: activation 1, inhibition = -1
			if ( tolower(type)=="activation" ) {
				weights = 1
			} else {
				weights = -1
			}
			# add all edges from controllers to controlleds.
			for(i_controllers in controllers) {
				for(i_controlleds in controlleds) {
					#skip edges to self
					if(i_controllers != i_controlleds) {
						#if edge already exists throw a warning and skip
						if(!(i_controlleds %in% unlist(graph::edges(mygraph, which=i_controllers)))) {
							mygraph = graph::addEdge(i_controllers, i_controlleds,mygraph,weights=weights)							
						} else {
							warning( paste("Problem while adding edge (weight=",weights,") to graph: Edge already exists between ",i_controllers," and ",i_controlleds," (weight=",unlist(graph::edgeWeights(mygraph))[paste(i_controllers,i_controlleds,sep=".")],"). Skipping this edge.",sep="" )  )
						}
					}
				}
			}			
			
		}
	}
	mygraph
}



#TODO transitive reduction & closure
#' This function generates the transitive closure of the supplied graph.
#' 
#' This function generates the transitive closure of the supplied graph. In short: if A->B->C then an edge A->C is added.
#' Edge weights are conserved if possible (in a hopefully smart way). 
#' 
#' @param mygraph graphNEL
#' @return Returns the transitive closure of the supplied graph.
#' @author Frank Kramer
#' @export
transitiveClosure <- function(mygraph) {
	
	if(!require(graph)) {
		message(paste("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!","\n"))
		return(NULL)
	}
	
	accessibles = acc(mygraph,graph::nodes(mygraph))
	for(n in names(accessibles)) {
		for(c in names(accessibles[[n]])) {
			if(accessibles[[n]][[c]]>1) {
				#cat(paste("adding",n,c,"\n"))
				mygraph = graph::addEdge(from=n,to=c, graph=mygraph, weights=1) ###FIX! XXX
			}
		}
	}
	mygraph
}

#' This function generates the transitive reduction of the supplied graph.
#' 
#' This function generates the transitive reduction of the supplied graph. In short: if A->B->C AND A->C then edge A->C is removed.
#' Edge weights are conserved if possible (in a hopefully smart way). 
#' 
#' @param mygraph graphNEL
#' @return Returns the transitive reduction of the supplied graph.
#' @author Frank Kramer
#' @export
transitiveReduction <- function(mygraph) {

	if(!require(graph)) {
		message(paste("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!","\n"))
		return(NULL)
	}
	
	temp=mygraph
	#remove temps wedge weights
	for(e in names(graph::edges(temp))) {
		for(t in graph::edges(temp)[[e]]) {
			graph::edgeData(temp,e,t,attr="weight") <-  1
		}
	}
	
	#enumarte pathws and remove all paths that have longer equivalents and same end node
	
#	ttt = as(mygraph,"graphAM")
#	as(ttt,"matrix")
#	# order by outdegree
#	# check crossing between every 2 edges?
#	
#	accessibles = acc(mygraph,nodes(mygraph))
#	for(n in names(accessibles)) {
#		for(c in names(accessibles[[n]])) {
#			if(accessibles[[n]][[c]]>1) {
#				graph::removeEdge(from=n,to=c, mygraph)
#			}
#		}
#	}
	mygraph
	
}

#' This function generates a (more or less) beautiful layout for a regulatory graph.
#' 
#' This function generates a (more or less) beautiful layout for a regulatory graph. 
#' Call this after you generated a graph with pathway2RegulatoryGraph. Since beauty is always in the eye of the beholder consider this a starting point for making your graphs even nicer.
#' Rgraphviz with dot layout is used.
#' Edges are green/red with normal/tee arrowheads for activations/inhibitions. 
#' If you want to specifically paint subgraphs in different colors use lists of vectors with node names for parameter subgraphs and vector of color names for subgraphs.color for your choice of color.   
#' The output can be further tweaked by setting layout options using nodeRenderInfo(mygraph) <- list() ...
#' See the Rgraphviz and Graphviz documentations.
#' 
#' @param mygraph graphNEL 
#' @param label Label of the graph
#' @param node.fixedsize logical. If font size is fixed or variable in regards to the nodes.
#' @param edge.weights vector. which colors to use for weighted edges
#' @param edge.arrowheads vector. which arrowheads to use for weighted edges
#' @param subgraphs A list of character vectors with node names defining the sub graphs.
#' @param subgraphs.colors vector. which colors to use for subgraphs
#' @return Returns the supplied graph in a layouted form with several parameters set for regulatory graph plotting.
#' @author Frank Kramer
#' @export
layoutRegulatoryGraph <- function(mygraph, label="", node.fixedsize=FALSE, edge.weights=c("green","black","red"), edge.arrowheads=c("normal","tee"),
		subgraphs=list(), subgraphs.colors=c("#B3E2CD","#FDCDAC","#F4CAE4","#E6F5C9","#FFF2AE")) {

	
	if(!require(graph)) {
		message("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!")
		return(NULL)
	}
	if(!require(Rgraphviz)) {
		message("This functions needs the Rgraphviz library installed, albeit it cannot be found. Check out the installation instructions!")
		return(NULL)
	}
	
	if(label != "")	graph::graphRenderInfo(mygraph) <- list(label=label, labelJust="l", labelLoc="t")
	graph::nodeRenderInfo(mygraph) <- list(shape="ellipse", fill="#e0e0e0", fixedsize=node.fixedsize)
	
	#graph::nodeRenderInfo(mygraph) <- list(iwidth="2", iheight="2")
	
	mygraph = Rgraphviz::layoutGraph(mygraph)		
	
	#GRAPH
	# graph rendering info
	if(label != "")	graph::graphRenderInfo(mygraph) <- list(label=label, labelJust="l", labelLoc="t")
	
	#EDGES
	# generate a list of all edges, preinitialized with the weights
	x = unlist(graph::edgeData(mygraph))
	names(x) = gsub(".weight","",names(x))
	names(x) = gsub("|","~",names(x), fixed=TRUE)
	
	#set edge arrowheads
	arrowhead = x
	arrowhead[TRUE] = rep(edge.arrowheads[1],length(x))
	arrowhead[x==1] = edge.arrowheads[1]
	arrowhead[x==-1] = edge.arrowheads[2]
	graph::edgeRenderInfo(mygraph) <- list(arrowhead=arrowhead)
	
	#set edge colors
	color = x
	color[TRUE] = rep(edge.weights[2],length(x))
	color[x==1] = edge.weights[1]
	color[x==-1] = edge.weights[3]
	graph::edgeRenderInfo(mygraph) <- list(col=color)
	
	#NODES
	#shape size and fill of nodes
	graph::nodeRenderInfo(mygraph) <- list(shape="ellipse", fill="#e0e0e0", fixedsize=node.fixedsize)
	
	#SUBGRAPHS
	if(length(subgraphs) > 0) {
		for(i in 1:length(subgraphs)) {
			for(n in subgraphs[[i]]) {
				graph::nodeRenderInfo(mygraph)$fill[n] <- subgraphs.colors[i]
			}
		}
	}
	
	return(mygraph)						
}

#' This function layouts a regulatory graph and plots it using Rgraphviz.
#' 
#' This function takes a regulatory graph as generated by pathway2regulatoryGraph and plots it using standard layout options of layoutRegulatoryGraph.
#' This function is a wrapper for layoutRegulatoryGraph with standard parameters.
#' Subgraphs can be painted with different colors. This can be done by passing parameter subgraph a list of character vectors with node names.
#' 
#' @param mygraph graphNEL, regulatory graph
#' @param subgraphs list of character vectors with node names
#' @param layoutGraph logical. If FALSE the graph is not laid out again but send directly to Rgraphviz::renderGraph.
#' @return none
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  pwid1 = "pid_p_100002_wntpathway"
#'  pwid2 = "pid_p_100146_hespathway"
#'  mygraph = pathway2RegulatoryGraph(biopax, pwid1)
#'  plotRegulatoryGraph(mygraph)
plotRegulatoryGraph <- function(mygraph, subgraphs=list(), layoutGraph=TRUE) {
	
	if(!require(Rgraphviz)) {
		message("This functions needs the Rgraphviz library installed, albeit it cannot be found. Check out the installation instructions!")
		return(NULL)
	}
	temp = mygraph
	if(layoutGraph) {
		temp = layoutRegulatoryGraph(temp, subgraphs=subgraphs)
	}
	Rgraphviz::renderGraph(temp)
}

#' This function calculates the overlap of 2 graphs
#' 
#' This function calculates the overlap of supplied graph1 with graph2.
#' Layout and weights of graph1 are kept.
#' 
#' @param graph1 graphNEL
#' @param graph2 graphNEL
#' @return Returns a list containing the compared graphs and edge- and node-wise overlap between them.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  pwid1 = "pid_p_100002_wntpathway"
#'  pwid2 = "pid_p_100146_hespathway"
#'  mygraph1 = pathway2RegulatoryGraph(biopax, pwid1)
#'  mygraph2 = pathway2RegulatoryGraph(biopax, pwid2)
#'  calcGraphOverlap(mygraph1,mygraph2)
calcGraphOverlap <- function(graph1, graph2) {

	if(!require(graph)) {
		message("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!")
		return(NULL)
	}
	
	edges_overlap = (length(graph::edgeNames(graph1)) - length(setdiff(graph::edgeNames(graph1), graph::edgeNames(graph2)))) / length(graph::edgeNames(graph1))
	nodes_overlap = (length(graph::nodes(graph1)) - length(setdiff(graph::nodes(graph1), graph::nodes(graph2)))) / length(graph::nodes(graph1))
	list(graph1=graph1, graph2=graph2, edges_overlap = edges_overlap , nodes_overlap = nodes_overlap)	
}

#' This function returns a graph computed by the insection of supplied graph1 and graph2.
#' 
#' This function returns a graph computed by the insection of supplied graph1 and graph2.
#' Layout and weights of graph1 are kept.
#' 
#' @param graph1 graphNEL
#' @param graph2 graphNEL
#' @return Returns the intersection of graph1 and graph2.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  pwid1 = "pid_p_100002_wntpathway"
#'  pwid2 = "pid_p_100146_hespathway"
#'  mygraph1 = pathway2RegulatoryGraph(biopax, pwid1)
#'  mygraph2 = pathway2RegulatoryGraph(biopax, pwid2)
#'  plotRegulatoryGraph(intersectGraphs(mygraph1,mygraph2))
intersectGraphs <- function(graph1, graph2) {
	
	if(!require(graph)) {
		message("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!")
		return(NULL)
	}
	
	temp = graph1
	e1 = strsplit(graph::edgeNames(graph1), split="~")
	e2 = strsplit(graph::edgeNames(graph2), split="~")
	i  = setdiff(e1,e2)
	for(x in 1:length(i)) {
		temp = graph::removeEdge(i[[x]][1],i[[x]][2], temp)
	}
	for(x in setdiff(graph::nodes(graph1), graph::nodes(graph2))) {
		temp = graph::removeNode(x, temp)
	}
	temp	
}

#' This function returns the different nodes and edges between graph1 and graph2.
#' 
#' This function returns the different nodes and edges between graph1 and graph2.
#' Layout options of graph1 are kept.
#' Coloring currently not implemented.
#' 
#' @param graph1 graphNEL
#' @param graph2 graphNEL
#' @param colorNodes logical
#' @param colors character vector of colors. If colorNodes==TRUE these colors are used for graph1 and graph2 respectivley.
#' @return Return the diff between the graphs.
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  pwid1 = "pid_p_100002_wntpathway"
#'  pwid2 = "pid_p_100146_hespathway"
#'  mygraph1 = pathway2RegulatoryGraph(biopax, pwid1)
#'  mygraph2 = pathway2RegulatoryGraph(biopax, pwid2)
#'  plotRegulatoryGraph(diffGraphs(mygraph1,mygraph2))
diffGraphs <- function(graph1, graph2, colorNodes=TRUE, colors=c("#B3E2CD","#FDCDAC") ) {
	
	if(!require(graph)) {
		message("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!")
		return(NULL)
	}
	
	temp = graph1
	e1 = strsplit(graph::edgeNames(graph1), split="~")
	e2 = strsplit(graph::edgeNames(graph2), split="~")
	i  = intersect(e1,e2)
	for(x in 1:length(i)) {
		temp = graph::removeEdge(i[[x]][1],i[[x]][2], temp)
	}
	temp
}

#' This function unites two graphs.
#' 
#' This function unites the two supplied graphs. Layout parameters from graph1 are used. 
#' If colorNodes==TRUE the returned graph has different colors for overlapping nodes and nodes individual for each graph. 
#' 
#' @param graph1 graphNEL
#' @param graph2 graphNEL
#' @param colorNodes logical
#' @param colors colors character vector of colors. If colorNodes==TRUE these colors are used for graph1 and graph2 respectivley.
#' @return Return a graph generated by uniting the two supplied graphs
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data
#'  data(biopax2example)
#'  pwid1 = "pid_p_100002_wntpathway"
#'  pwid2 = "pid_p_100146_hespathway"
#'  mygraph1 = pathway2RegulatoryGraph(biopax, pwid1)
#'  mygraph2 = pathway2RegulatoryGraph(biopax, pwid2)
#'  plotRegulatoryGraph(uniteGraphs(mygraph1,mygraph2))
uniteGraphs <- function(graph1, graph2, colorNodes=TRUE, colors=c("#B3E2CD","#FDCDAC","#F4CAE4")) {   #color: yellow, lightblue, lightgreen
	
	if(!require(graph)) {
		message("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!")
		return(NULL)
	}
	if(!require(Rgraphviz)) {
		message("This functions needs the Rgraphviz library installed, albeit it cannot be found. Check out the installation instructions!")
		return(NULL)
	}
	
	temp = graph1
	graph::graphRenderInfo(temp)$laidout <- FALSE
	
	#ADD NODES
	for(i in setdiff(graph::nodes(graph2), graph::nodes(graph1))) {
		temp = graph::addNode(i, temp)
	}
	
	#ADD EDGES
	e1 = strsplit(graph::edgeNames(graph1), split="~")
	e2 = strsplit(graph::edgeNames(graph2), split="~")
	i  = setdiff(e2,e1)
	for(x in 1:length(i)) {
		temp = graph::addEdge(i[[x]][1],i[[x]][2], temp, as.numeric(graph::edgeData(graph2,i[[x]][1],i[[x]][2], "weight")))
	}

	#shape size and fill of nodes
	nodevector = as.vector(rep(1,length(graph::nodes(temp))))
	names(nodevector) = graph::nodes(temp)
	
	node.shape = nodevector
	node.shape[TRUE] = "ellipse"
	graph::nodeRenderInfo(temp)$shape <- node.shape
	
	#LAYOUT GRAPH
	temp = Rgraphviz::layoutGraph(temp)
	
	#set shape again. strange rgraphviz behavior
	graph::nodeRenderInfo(temp)$shape <- node.shape

	node.fill = nodevector
	node.fill[TRUE] = colors[1]
	graph::nodeRenderInfo(temp)$fill <- node.fill

	### reuse previous layout of graph1
	if(!is.null(graph::nodeRenderInfo(graph1)$height[1])) {
		node.height = nodevector
		node.height[TRUE] = graph::nodeRenderInfo(graph1)$height[1]
		graph::nodeRenderInfo(temp)$height <- node.height	
	}
	if(!is.null(graph::nodeRenderInfo(graph1)$width[1])) {
		node.width = nodevector
		node.width[TRUE] = graph::nodeRenderInfo(graph1)$width[1]
		graph::nodeRenderInfo(temp)$width <- node.width
	}
	if(!is.null(graph::nodeRenderInfo(graph1)$fontsize[1])) {
		node.fontsize = nodevector
		node.fontsize[TRUE] = graph::nodeRenderInfo(graph1)$fontsize[1]
		graph::nodeRenderInfo(temp)$fontsize <- node.fontsize
	}
	if(!is.null(graph::nodeRenderInfo(graph1)$labelfontsize[1])) {
		node.labelfontsize = nodevector
		node.labelfontsize[TRUE] = graph::nodeRenderInfo(graph1)$labelfontsize[1]
		graph::nodeRenderInfo(temp)$labelfontsize <- node.labelfontsize
	}
	if(!is.null(graph::nodeRenderInfo(graph1)$labelfontsize[1])) {
		node.fixedsize = nodevector
		node.fixedsize[TRUE] = graph::nodeRenderInfo(graph1)$fixedsize[1]
		graph::nodeRenderInfo(temp)$fixedsize <- node.fixedsize
	}
	
	#COLOR SUBGRAPHS
	#only graph1
	for(n in setdiff(graph::nodes(graph1), graph::nodes(graph2))) {
		graph::nodeRenderInfo(temp)$fill[n] <- colors[2]
	}
	#only graph2
	for(n in setdiff(graph::nodes(graph2), graph::nodes(graph1))) {
		graph::nodeRenderInfo(temp)$fill[n] <- colors[3]
	}
	
	#COLOR EDGES	
	# generate a list of all edges, preinitialized with the weights
	edge.weights=c("green","black","red")
	x = unlist(graph::edgeData(temp))
	names(x) = gsub(".weight","",names(x))
	names(x) = gsub("|","~",names(x), fixed=TRUE)
	
	#set edge arrowheads
	edge.arrowheads=c("normal","tee")
	arrowhead = x
	arrowhead[TRUE] = rep(edge.arrowheads[1],length(x))
	arrowhead[x==1] = edge.arrowheads[1]
	arrowhead[x==-1] = edge.arrowheads[2]
	graph::edgeRenderInfo(temp) <- list(arrowhead=arrowhead)
	
	#set edge colors
	color = x
	color[TRUE] = rep(edge.weights[2],length(x))
	color[x==1] = edge.weights[1]
	color[x==-1] = edge.weights[3]
	graph::edgeRenderInfo(temp) <- list(col=color)
	
	return(temp)
}

#' This function colors the nodes of a graph.
#' 
#' This function colors nodes of a graph, usually this is used to color subgraphs 
#' or add a color hue correlating with the expression level or fold change to the molecules. 
#' 
#' @param graph1 graphNEL
#' @param nodes vector of node names specifiying which nodes to color. must be same length as parameter foldChanges
#' @param values vector of values indicating fold changes, gene expression values or similar. colors are mapped linearly over the range of these values
#' @param colors string. either "greenred" or "yellowred", specifying which color gradient to use.
#' @return Returns a graph with specified nodes colored according to the foldChanges
#' @author Frank Kramer
#' @export
#' @examples
#'  # load data and retrieve wnt pathway
#'  data(biopax2example)
#'  pwid1 = "pid_p_100002_wntpathway"
#'  mygraph1 = pathway2RegulatoryGraph(biopax, pwid1)
#'  mygraph1 = layoutRegulatoryGraph(mygraph1)
#'  # retrieve all nodes
#'  nodes = nodes(mygraph1)
#'  # random expression data for your nodes 
#'  values = rnorm(length(nodes), mean=6, sd=2)
#'  # color nodes of the graph
#'  mygraph1 = colorGraphNodes(mygraph1, nodes, values, colors="greenred") 
#'  # plot the now colored graph 
#'  plotRegulatoryGraph(mygraph1, layoutGraph=FALSE)
colorGraphNodes <- function(graph1, nodes, values, colors=c("greenred","yellowred")) {
	
	if(!require(graph)) {
		message(paste("This functions needs the graph library installed, albeit it cannot be found. Check out the installation instructions!","\n"))
		return(NULL)
	}
	
	ncol = length(nodes)
	r=seq(0,1,length=ncol)
	r2=seq(-1,1,length=ncol)
	if(colors=="greenred") {
		col = hcl(h=ifelse(r2>0, 90-90*abs(r2), 90+90*abs(r2)), c=30+50*abs(r2)^0.2, l=92-72*abs(r2)^1.5)
	}
	if(colors=="yellowred") {
		col = hcl(h=90-90*r, c=30+50*r^0.2, l=90-60*r^1.5)
	}
	lim = range(values)
	x=findInterval(values, seq(lim[1],lim[2],diff(lim)/ncol)[-c(1,ncol+1)])+1
	res=col[x]
	
	for(n in 1:length(nodes)) {
		graph::nodeRenderInfo(graph1)$fill[nodes[n]] = res[n] 
	}
	 
	return(graph1)
}
	
	
	
#	
#	
#	
#	
#	
#	########## OLD STUFF ################
#	graph::graphRenderInfo(temp)$laidout <- FALSE
#	
#	
#	#shape size and fill of nodes
#	nodevector = as.vector(rep(1,length(graph::nodes(graph1))))
#	names(nodevector) = graph::nodes(graph1)
#	
#	node.fill = nodevector
#	if(nodes == NA)	{
#		node.fill[TRUE] = colors
#		graph::nodeRenderInfo(graph)$fill <- node.fill
#	} else {
#		
#	}
#	
#	graph::nodeRenderInfo(temp)$shape <- node.shape
#	graph::nodeRenderInfo(temp)$fill <- node.fill
#	graph::nodeRenderInfo(temp)$height <- node.height
#	graph::nodeRenderInfo(temp)$width <- node.width
#	graph::nodeRenderInfo(temp)$fontsize <- node.fontsize
#	graph::nodeRenderInfo(temp)$labelfontsize <- node.labelfontsize
#	graph::nodeRenderInfo(temp)$fixedsize <- node.fixedsize
#	
#	#COLOR SUBGRAPHS
#	#only graph1
#	for(n in setdiff(graph::nodes(graph1), graph::nodes(graph2))) {
#		graph::nodeRenderInfo(temp)$fill[n] <- colors[2]
#	}
#	#only graph2
#	for(n in setdiff(graph::nodes(graph2), graph::nodes(graph1))) {
#		graph::nodeRenderInfo(temp)$fill[n] <- colors[3]
#	}
#	
#	#EDGES	
#	e1 = strsplit(graph::edgeNames(graph1), split="~")
#	e2 = strsplit(graph::edgeNames(graph2), split="~")
#	i  = setdiff(e2,e1)
#	for(x in 1:length(i)) {
#		temp = graph::addEdge(i[[x]][1],i[[x]][2], temp, as.numeric(graph::edgeData(graph2,i[[x]][1],i[[x]][2], "weight")))
#	}
#	# generate a list of all edges, preinitialized with the weights
#	edge.weights=c("green","black","red")
#	x = unlist(graph::edgeData(temp))
#	names(x) = gsub(".weight","",names(x))
#	names(x) = gsub("|","~",names(x), fixed=TRUE)
#	
#	#set edge arrowheads
#	edge.arrowheads=c("normal","tee")
#	arrowhead = x
#	arrowhead[TRUE] = rep(edge.arrowheads[1],length(x))
#	arrowhead[x==1] = edge.arrowheads[1]
#	arrowhead[x==-1] = edge.arrowheads[2]
#	graph::edgeRenderInfo(temp) <- list(arrowhead=arrowhead)
#	
#	#set edge colors
#	color = x
#	color[TRUE] = rep(edge.weights[2],length(x))
#	color[x==1] = edge.weights[1]
#	color[x==-1] = edge.weights[3]
#	graph::edgeRenderInfo(temp) <- list(col=color)
#	
#	Rgraphviz::layoutGraph(temp)	
#}


