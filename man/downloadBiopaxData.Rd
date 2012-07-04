\name{downloadBiopaxData}
\alias{downloadBiopaxData}
\title{This function downloads Biopax data from online databases}
\usage{
  downloadBiopaxData(database = "NCI",
    model = c("pid", "biocarta", "reactome", outputfile = ""),
    version = "biopax2")
}
\arguments{
  \item{database}{string. Select which database you want to
  download from. Currently only NCI links have been
  stored.}

  \item{model}{string. Select which model/file you want to
  download. Currently NCI versions of the Pathway
  Interaction Database, Biocarta and Reactome are linked.}

  \item{version}{string. Select which Biopax Version you
  want to download.}

  \item{outputfile}{string. The file name to save the
  downloaded data in. If left empty the URL file name will
  be used. The unzipped file name can be different from
  this. Check the screen output of gunzip.}
}
\value{
  none. Check output for the name of the unzipped biopax
  .owl file.
}
\description{
  This function has an internal list of download links for
  some online databases. It will retrieve the selected
  model from the selected database using RCurl. The
  downloaded file is (if needed) unzipped and ready to be
  used as input for rBiopaxParser::readBiopax.
}
\author{
  fkramer
}

