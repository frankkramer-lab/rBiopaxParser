###############################################################################
#
# test_downloadBiopaxData.R: 	This file contains all RUnit tests for functions in downloadBiopaxData.R.
# author: Frank Kramer <dev@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################

test_DATABASE_BIOPAX <- function() {
	checkTrue( "data.frame" == class(DATABASE_BIOPAX) )
} 

test_downloadBiopaxData <- function() {
	return(TRUE)
} 
