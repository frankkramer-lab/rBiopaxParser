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
	checkEquals(getSubClasses("bp:control"), getSubClasses("control"))
	
	#check that "control" is a subclass of "interaction" in BP2
	checkTrue( "control" %in% getSubClasses("interaction") )
	
	#check that "dna" is a subclass of "physicalEntity" in BP2
	checkTrue( "dna" %in% getSubClasses("physicalEntity") )
	
} 

test_getSuperClasses <- function() {
	
	#check that namespaces are correctly stripped off class names
	checkEquals(getSuperClasses("bp:control"), getSuperClasses("control"))
	
	#check that "interaction" is a superclass of "control" in BP2
	checkTrue( "interaction" %in% getSuperClasses("control") )
	
	#check that "physicalEntity" is a superclass of "dna" in BP2
	checkTrue( "physicalEntity" %in% getSuperClasses("dna") )
	
} 

test_getClassProperties <- function() {
	
	#check that namespaces are correctly stripped off class names
	checkEquals(getClassProperties("bp:control"), getClassProperties("control"))
	
	#check that "control" has a property "CONTROL-TYPE" in BP2
	checkTrue( "CONTROL-TYPE" %in% unlist(getClassProperties("control")$property) )
	
} 

