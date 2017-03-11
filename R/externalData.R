#' SNV data.
#'
#' A MAF formattted file containing TCAG OV carcinoma SNVs for 16
#' genes.
#'
#' @format A MAF formatted file with 15 rows and 32 columns:
#' \describe{
#'   \item{HUGO_Symbol}{HGNC gene symbol}
#'   \item{ENTREZ_GENE_ID}{ID}
#'   ...
#' }
#' @source \url{https://tcga-data.nci.nih.gov/docs/publications/ov_2011/}
#' @examples
#' \dontrun{
#' system.file("extdata", "devSNV.maf", package="rete")
#' fPath <- system.file("extdata", "devSNV.maf", package="rete")
#' myMAF <- read.delim(fPath, skip = 1, stringsAsFactors = FALSE)
#' }
#' @docType data
#' @name devSNV.maf
NULL

#' CNA data.
#'
#' A file containing TCAG OV carcinoma CNAs for 16
#' genes.
#'
#' @format A tab separated text file with 16 rows and 582 columns, most of
#'   which represent individual cancer tissue samples:
#' \describe{
#'   \item{Gene Symbol}{HGNC gene symbol}
#'   \item{Locus ID}{Entrez Gene ID}
#'   \item{Cytoband}{Chromosome, arm and band information}
#'   \item{<Sample Barcode>...}{497 tissue sample barcodes}
#'   ...
#' }
#' @source \url{https://tcga-data.nci.nih.gov/docs/publications/ov_2011/}
#' @examples
#' \dontrun{
#' system.file("extdata", "devCNA.txt", package="rete")
#' fPath <- system.file("extdata", "devCNA.txt", package="rete")
#' myCNA <- read.delim(fPath, stringsAsFactors = FALSE)
#' }
#' @docType data
#' @name devCNA.txt
NULL

#' PPI data.
#'
#' A file containing protein-protein interactions for 16
#' genes, subset from iRefIndex version 14 data.
#'
#' @format A MITAB formatted text file with 75 rows and 54 columns:
#' \describe{
#'   \item{uidA}{UniProt ID of interactor A}
#'   \item{uidB}{UniProt ID of interactor B}
#'   ...
#' }
#' @source \url{http://irefindex.org/}
#' @examples
#' \dontrun{
#' system.file("extdata", "devPPI.txt", package="rete")
#' fPath <- system.file("extdata", "devPPI.txt", package="rete")
#' myMITAB <- read.delim(fPath, stringsAsFactors = FALSE)
#' }
#' @docType data
#' @name devPPI.txt
NULL

