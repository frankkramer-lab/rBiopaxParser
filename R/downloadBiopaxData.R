###############################################################################
#
# downloadBiopaxData.R: 	This file contains all functions related to retrieving Biopax data from the web.
# author: Frank Kramer <dev@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################

#' Databases available for direct download via downloadBiopaxData
#' 
#' A data.frame listing all available databases which can be directly downloaded (Homo Sapiens only) via function downloadBiopaxData.
#' The variables are as follows:
#' 
#' \itemize{
#'   \item database. Name of the database
#'   \item model. Name of the ontology model
#'   \item version. Biopax level
#'   \item link. Link to the direct download
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name DATABASE_BIOPAX
#' @title DATABASE_BIOPAX
#' @usage DATABASE_BIOPAX
#' @format A data frame with 46 rows and 4 columns
#' @export
DATABASE_BIOPAX = data.frame(
		matrix(ncol=4,byrow=T, dimnames=list(list(),list("database","model","version","link")),data= list(
						"NCI",			"pid",		"biopax2",		"ftp://ftp1.nci.nih.gov/pub/PID/BioPAX_Level_2/NCI-Nature_Curated.bp2.owl.gz",					
						"NCI",			"biocarta",	"biopax2",		"ftp://ftp1.nci.nih.gov/pub/PID/BioPAX_Level_2/NCI-Nature_Curated.bp2.owl.gz",
						"NCI",			"reactome",	"biopax2",		"ftp://ftp1.nci.nih.gov/pub/PID/BioPAX_Level_2/NCI-Nature_Curated.bp2.owl.gz",
						"NCI",			"kegg",		"biopax2",		"ftp://ftp1.nci.nih.gov/pub/PID/BioPAX_Level_2/NCI-Nature_Curated.bp2.owl.gz",
						
						"NCI",			"pid",		"biopax3",		"http://sourceforge.net/mailarchive/attachment.php?list_name=biopax-paxtools&message_id=CALXvGpbptVgL2YRfK6VCwQpHH7hx3tOqZgnSYU9Uiff-8OWSLw%40mail.gmail.com&counter=3",					
						"NCI",			"biocarta",	"biopax3",		"http://sourceforge.net/mailarchive/attachment.php?list_name=biopax-paxtools&message_id=CALXvGpbptVgL2YRfK6VCwQpHH7hx3tOqZgnSYU9Uiff-8OWSLw%40mail.gmail.com&counter=4",
						"NCI",			"reactome",	"biopax3",		"http://sourceforge.net/mailarchive/attachment.php?list_name=biopax-paxtools&message_id=CALXvGpbqbVtj30uw5D2LmAgpV_fdLbDzDNCShBuwa-a%2BwTRy%3Dw%40mail.gmail.com&counter=2",
						"NCI",			"kegg",		"biopax3",		"ftp://ftp1.nci.nih.gov/pub/PID/BioPAX_Level_3/KEGG.bp3.owl.gz",
						
						"reactome",		"reactome",	"biopax2",		"http://www.reactome.org/download/current/biopax2.zip",
						"reactome",		"reactome",	"biopax3",		"http://www.reactome.org/download/current/biopax.zip"
				)), stringsAsFactors=FALSE)


#' This function downloads Biopax data from online databases
#' 
#' This function has an internal list of download links for some online databases. It will retrieve the selected model from the selected database using RCurl.
#' The downloaded file is (if needed) unzipped and ready to be used as input for rBiopaxParser::readBiopax.
#' This function requires package RCurl to run. 
#' You can easily skip this step by downloading the exported file yourself and continuing with readBiopax. 
#' 
#' @param database string. Select which database you want to download from. Currently only NCI links have been stored.
#' @param model string. Select which model/file you want to download. Currently NCI versions of the Pathway Interaction Database, Biocarta, Reactome and KEGG are linked.
#' @param version string. Select which Biopax Version you want to download.
#' @param outputfile string. The file name to save the downloaded data in. If left empty the URL file name will be used. The unzipped file name can be different from this. Check the screen output of gunzip.  
#' @return none. Check output for the name of the unzipped biopax .owl file.
#' @author fkramer
#' @export
#' @examples
#'  \dontrun{file = downloadBiopaxData("NCI", "biocarta", version = "biopax2")}
#'  \dontrun{biopax = readBiopax(file)}
#'  \dontrun{biopax}
downloadBiopaxData <- function(database="NCI", model=c("pid","biocarta","reactome", "kegg"), outputfile="", version="biopax2") {
	
	links = DATABASE_BIOPAX
		
	link = links[links$database==database & links$model==model[1] & links$version==version,"link"][[1]]
	
	m <- regexec("^(([^:]+)://)?([^:/]+)(:([0-9]+))?(/.*/)?(.*)", link)
	filename = regmatches(link, m)[[1]][length(m[[1]])]

	if(!require(RCurl)) {
		message(paste("This functions needs the RCurl library installed, albeit it cannot be found. Check out the installation instructions or manually download the file you wanted at: ",filename))
		return(0)
	}
	
	message(paste("Trying to download:",filename))
	if(database=="reactome") message("This download is very large (60-80MB) and might take a while.\n")
	content = RCurl::getBinaryURL(link)
	if(outputfile=="") outputfile = filename
	writeBin(content, useBytes = TRUE, con = outputfile )
	
	# for reactome: extract homo sapiens data
	if(database=="reactome") {
		content2 = readBin(unz(outputfile, "Homo sapiens.owl", open="rb"),"raw", n=500000000)
		outputfile = paste(outputfile,".owl",sep="")
		writeBin(content2, useBytes = TRUE, con = outputfile )
	}
	
	message(paste("Downloaded ", outputfile, "! Proceed with readBiopax!\n", sep=""))
	outputfile
}


# x = downloadBiopaxData(model="biocarta")