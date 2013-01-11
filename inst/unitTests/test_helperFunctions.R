###############################################################################
#
# test_helperFunctions.R:  This file contains all RUnit tests for functions in helperFunctions.R.
# author: Frank Kramer <dev@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################


test_unfactorize <- function() {
	
	#df1 with factor levels
	df1 = data.frame(num=c(1,2,3,4,5), string=c("one","two","three","four","five"), stringsAsFactors=TRUE)
	#df2 without factors
	df2 = data.frame(num=c(1,2,3,4,5), string=c("one","two","three","four","five"), stringsAsFactors=FALSE)
	
	#check that factors are correctly stripped off data.frames
	checkEquals(rBiopaxParser:::unfactorize(df1), df2)

} 

test_isOfNamespace <- function() {
	
	#check that namespace is correctly identified
	checkTrue( rBiopaxParser:::isOfNamespace(classname="bp:control",namespace="bp") )
	
} 

test_addns <- function() {
	
	#check that namespace is correctly added
	checkEquals( "bp:control", rBiopaxParser:::addns("control","bp") )
	
} 

test_stripns <- function() {
	
	#check that namespace is correctly stripped
	checkEquals( "control", rBiopaxParser:::stripns("bp:control") )
	
} 

test_addhash <- function() {
	
	#check that hash is correctly added
	checkEquals( "#id1", rBiopaxParser:::addhash("id1") )
	
} 

test_striphash <- function() {
	
	#check that hash is correctly stripped
	checkEquals( "id1", rBiopaxParser:::striphash("#id1") )
	
} 

