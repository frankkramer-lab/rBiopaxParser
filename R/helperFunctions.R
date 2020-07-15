###############################################################################
#
# helperFunctions.R: 	This file contains the all helper functions not directly related to any other source file.
# author: Frank Kramer <dev@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################

#' Replace factors/levels in a data.frame and use plain strings instead 
#' 
#' This function takes a data.frame as argument and returns it with strings instead of factors.
#' 
#' @param df any data.frame with factor levels in at least one column
#' @return The data.frame is returned using strings instead of factors.
#' @author Frank Kramer
unfactorize <- function (df) {
	if(!("data.frame" %in% class(df))) { stop("Error: unfactorize: data.frame as argument expected") }
	
	for (i in 1:ncol(df)) {
		if (class(df[,i]) == "factor") df[,i] <- as.character(df[,i])
	}
	df
}


################### HELPER FUNCTIONS START ########
#class and property checking

#' Check if a classname is preceeded by a certain namespace tag like in "namespace:classname"
#' 
#' This function checks if the supplied input string starts with a supplied namespace tag
#' 
#' @param classname A string containing the classname to check
#' @param namespace A string giving the namespace to check for
#' @return This function returns TRUE if the supplied classname string is preceeded with the supplied namespace string, and FALSE if not.
#' @author Frank Kramer
isOfNamespace <- function(classname, namespace="bp") {
	if(is.null(namespace) | is.na(namespace) | namespace == "" ) {
		!grepl(":",classname,ignore.case = FALSE)
	} else {
		grepl(paste(namespace,":",sep=""),classname,ignore.case = FALSE)
	}
} 

#' Add a namespace tag to the supplied classname string
#' 
#' This function takes the input classname, checks if it already has a namespace, and if not pastes the namespace tag with a dividing ":" in front of it.
#' 
#' @param classname A string containing a classname
#' @param namespace A string containing a namespace
#' @return If the classname is not preceeded by a namespace yet, the supplied namespace is pasted in front of it and returned.
#' @author Frank Kramer
addns <- function(classname, namespace="bp") {
	sel=!is.na(classname) & !(nchar(classname)==0) & !grepl(":",classname,ignore.case = FALSE)
	classname[sel] = paste(namespace,":",classname[sel],sep="")
	classname
}


#' Strips a namespace tag off a supplied classname string
#' 
#' Strips a namespace tag off a supplied classname string
#'  
#' @param classname A string containing a classname preceeded by a namespace tag
#' @return The classname with the namespace tag stripped off it.
#' @author Frank Kramer
stripns <- function(classname) {
	gsub(".*:", "", classname, perl=TRUE)
}

#' Check if a string is an URL, preceeded by "http:"
#' 
#' This function checks if the supplied input string starts with "http:"
#' 
#' @param string A string containing the classname to check
#' @return This function returns TRUE if the supplied classname string starts with "http:", and FALSE if not.
#' @author Frank Kramer
isURL <- function(string) {
	grepl(pattern="http:", string, ignore.case=T)
} 

#' Adds a hash in front of a string
#' 
#' Adds a hash in front of a string 
#' 
#' @param x A string to be preceeded by a hash
#' @return The supplied string with a hash "#" pasted in front of it.
#' @author Frank Kramer
addhash <- function(x) {
	paste("#",x,sep="")
}

#' Strips a hash in front of a string
#' 
#' Strips a hash in front of a string
#' @param x A string to be stripped off a preceeeding hash
#' @return The supplied string with a hash "#" stripped off front.
#' @author Frank Kramer
striphash <- function(x) {
	sub("#","",x)
}

#' Checks if instances in the biopax data.table are of the given class
#' 
#' This function checks if instances in the supplied biopax data.table are of a given class. 
#' If considerInheritance is set to TRUE it also checks if instances are of a given class or any of its inherited classes.
#' 
#' @param df A data.frame with biopax instances 
#' @param class A string containing the class name to check for
#' @param considerInheritance Logical value indicating wether to consider inheritance or not
#' @param biopaxlevel Numeric. Specifies the Biopax Level to use.
#' @return Returns TRUE for every row in the data.frame which is of the supplied class 
#' @author Frank Kramer
#' @import data.table
#' @export
#' @examples
#'  # load data
#'  data(biopaxexample)
isOfClass <- function (df, class, considerInheritance=FALSE, biopaxlevel=2) {
	class = stripns(class)
	if(considerInheritance) {
		class = c(class,getSubClasses(class,biopaxlevel))
	}
	tolower(as.character(df$class)) %chin% tolower(class)
}

#' Checks if instances in the biopax data.table have a given property
#' 
#' Checks if instances in the biopax data.table have a given property
#'  
#' @param df A data.frame with biopax instances 
#' @param property A string containing the name of the property to check for
#' @return Returns TRUE for every row in the data.frame with contains the supplied property. Logical vector with length corresponding to the number of rows in the data.frame.
#' @author Frank Kramer
#' @import data.table
#' @examples
#'  # load data
#'  data(biopaxexample)
#' @export
hasProperty <- function (df, property) {
	tolower(as.character(df$property)) %chin% stripns(tolower(property))
}

#' This function checks the supplied arguments if they abid to the given restrictions
#' 
#' This function checks the supplied arguments if they abid to the given restrictions 
#' 
#' @param args The vector of arguments to check
#' @param allowedValues A named list of values the argument of a this name is allowed to have
#' @param allowNULL Logical, allow NULL or not 
#' @param allowNA Logical, allow NA or not
#' @param allowEmptyString Logical, allow empty strings or not
#' @param allowInf Logical, allow values of +/- infinity or not
#' @return Returns 1 if all checks completed successfully, returns error message otherwise.
#' @author Frank Kramer
internal_checkArguments <- function(args = c(), allowedValues=list(), allowNULL=FALSE, allowNA=FALSE, allowEmptyString=TRUE, allowInf=TRUE  ) {
	
	# return value == 1 -> everything ok, otherwise error message is returned
	# first check if call to this function is ok
	if(class(allowNA) != "logical") stop("Function .checkArguments expects parameter 'allowNA' to be of class 'logical'")
	if(class(allowNULL) != "logical") stop("Function .checkArguments expects parameter 'allowNA' to be of class 'logical'")
	if(class(allowEmptyString) != "logical") stop("Function .checkArguments expects parameter 'allowNA' to be of class 'logical'")
	if(class(allowInf) != "logical") stop("Function .checkArguments expects parameter 'allowNA' to be of class 'logical'")
	if(!is.list(args)) stop("Function .checkArguments expects parameter 'args' to be of class 'list'")
	if(length(args) == 0) stop("Function .checkArguments expects parameter 'args' to have length > 0")
	if( any(sapply(args,class) != "character") ) stop("Function .checkArguments expects parameter 'args' to be a list of strings.")
	if(!is.list(allowedValues)) stop("Function .checkArguments expects parameter 'allowedValues' to be of class 'list'")
	if(length(allowedValues) != 0 & length(allowedValues) != length(args)) stop("Function .checkArguments expects parameter 'allowedValues' to have length = 0 or same length as parameter 'args'")
	
	args = unlist(args)
	#check passed arguments for existance, null, na, empty string, infinite, class
	for(i in 1:length(args)) {
		if(!exists(args[i])) { return(paste("Error: Variable ", args[i], " is undefined.",sep="")) }
		if(!allowNULL) { if( is.null(get(args[i])) ) { return(paste("Error: Variable ", args[i], " is null.",sep="")) }  }
		if(!allowNA) { if( is.na(get(args[i])) ) { return(paste("Error: Variable ", args[i], " is na.",sep="")) }  }
		if(!allowEmptyString) { if( get(args[i]) == "" ) { return(paste("Error: Variable ", args[i], " is an empty string.",sep="")) }  }
		if(!allowInf) { if( is.infinite(get(args[i])) ) { return(paste("Error: Variable ", args[i], " is infinite.",sep="")) }  }
		
		if(args[i] %in% names(allowedValues)) {
			if( !(any(class(get(args[i])) %in% allowedValues[args[i]])) ) {
				return(paste("Error: Variable ", args[i], " is not of class ",paste(allowedValues[args[i]], collapse=" "),".",sep=""))
			}
		}
	}
	return(1)
}


