#' Parses BioPax files and represents them in R
#'
#' rBiopaxParser is a...
#'
#' \tabular{ll}{
#' Package: \tab rBiopaxParser\cr
#' Type: \tab Package\cr
#' Version: \tab 0.15\cr
#' Date: \tab 2012-08-22\cr
#' License: \tab GPL (>= 2)\cr
#' }
#'
#' @author Frank Kramer \email{dev@@frankkramer.de}
#' @name rBiopaxParser-package
#' @aliases rBiopaxParser
#' @docType package
#' @title Parses BioPax level files and represents them in R
#' @keywords package
#' @examples
#' \dontrun{biopax = readBiopax(file="biopaxmodel.owl")}
NULL

#' Biopax example data set
#' 
#' A dataset containing two regulatory pathways encoded in Biopax Level 2 and parsed in via readBiopax().
#' 
#' @docType data
#' @keywords datasets
#' @format An example biopax model parsed in via readBiopax.
#' @name biopax
#' @alias biopax2example
#' @examples
#' data(biopax2example)
#' biopax
NULL