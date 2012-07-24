###############################################################################
#
# downloadBiopaxData.R: 	This file contains all functions related to retrieving Biopax data from the web.
# author: Frank Kramer <mail@frankkramer.de>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################

#' This function downloads Biopax data from online databases
#' 
#' This function has an internal list of download links for some online databases. It will retrieve the selected model from the selected database using RCurl.
#' The downloaded file is (if needed) unzipped and ready to be used as input for rBiopaxParser::readBiopax. 
#' 
#' @param database string. Select which database you want to download from. Currently only NCI links have been stored.
#' @param model string. Select which model/file you want to download. Currently NCI versions of the Pathway Interaction Database, Biocarta and Reactome are linked.
#' @param version string. Select which Biopax Version you want to download.
#' @param outputfile string. The file name to save the downloaded data in. If left empty the URL file name will be used. The unzipped file name can be different from this. Check the screen output of gunzip.  
#' @returnType none
#' @return none. Check output for the name of the unzipped biopax .owl file.
#' @author fkramer
#' @export
downloadBiopaxData <- function(database="NCI", model=c("pid","biocarta","reactome"), outputfile="", version="biopax2") {
	
	links = data.frame(database=c("NCI"), model=c("pid"), version=c("biopax2"), link=c("ftp://ftp1.nci.nih.gov/pub/PID/BioPAX_Level_2/NCI-Nature_Curated.bp2.owl.gz"), stringsAsFactors=FALSE)
	links = rbind(links,c("NCI","biocarta","biopax2","ftp://ftp1.nci.nih.gov/pub/PID/BioPAX_Level_2/BioCarta.bp2.owl.gz"))
	links = rbind(links,c("NCI","reactome","biopax2","ftp://ftp1.nci.nih.gov/pub/PID/BioPAX_Level_2/Reactome.bp2.owl.gz"))
	
	link = links[links$database==database & links$model==model[1] & links$version==version,"link"]
	
	m <- regexec("^(([^:]+)://)?([^:/]+)(:([0-9]+))?(/.*/)?(.*)", link)
	filename = regmatches(link, m)[[1]][length(m[[1]])]

	if(!require(RCurl)) {
		cat(paste("This functions needs the RCurl library installed, albeit it cannot be found. Check out the installation instructions or manually download the file you wanted at: ",filename,"\n"))
		return(0)
	}
	
	cat(paste("Trying to download:",filename,"\n"))
	content = RCurl::getBinaryURL(link)
	if(outputfile=="") outputfile = filename
	writeBin(content, useBytes = TRUE, con = outputfile )
	
	cat("Unzipping using gunzip \n")
	system(paste("gunzip --verbose",outputfile))
}


# x = downloadBiopaxData(model="biocarta")