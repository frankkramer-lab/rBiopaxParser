###############################################################################
#
# test_biopax2Classes.R: This file contains all RUnit tests for functions in biopax2Classes.R.
# author: Frank Kramer <dev@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################


test_DATABASE_BIOPAX <- function() {
	checkTrue( "data.frame" == class(CLASS_INHERITANCE_BP2) )
} 

test_DATABASE_BIOPAX <- function() {
	checkTrue( "data.frame" == class(CLASS_INHERITANCE_BP3) )
} 

test_DATABASE_BIOPAX <- function() {
	checkTrue( "data.frame" == class(CLASS_PROPERTIES_BP2) )
} 

test_DATABASE_BIOPAX <- function() {
	checkTrue( "data.frame" == class(CLASS_PROPERTIES_BP3) )
} 

test_getSubClasses <- function() {
	
	#check that namespaces are correctly stripped off class names
	checkEquals(getSubClasses("bp:Control"), getSubClasses("Control"))
	
	#check that "control" is a subclass of "interaction"
	checkTrue( "control" %in% tolower(getSubClasses("Interaction")) )
	
	#check that "dna" is a subclass of "PhysicalEntity"
	checkTrue( "dna" %in% tolower(getSubClasses("PhysicalEntity")) )
	
} 

test_getSuperClasses <- function() {
	
	#check that namespaces are correctly stripped off class names
	checkEquals(getSuperClasses("bp:Control"), getSuperClasses("Control"))
	
	#check that "interaction" is a superclass of "control"
	checkTrue( "interaction" %in% tolower(getSuperClasses("Control")) )
	
	#check that "physicalEntity" is a superclass of "dna"
	checkTrue( "physicalentity" %in% tolower(getSuperClasses("Dna")) )
	
} 

test_getClassProperties <- function() {
	
	#check that namespaces are correctly stripped off class names
	checkEquals(getClassProperties("bp:Control"), getClassProperties("Control"))
	
	#check that "control" has a property "CONTROL-TYPE"
	checkTrue( tolower("controlType") %in% tolower(getClassProperties("Control")$property) )
	
} 

