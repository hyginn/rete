# importCNA.GISTIC2.R
#
#
# Imports copy number aberration data from Firehose database,
# produces an rCNA object saved as a RDS file.
#
#
#'@param fFHG A vector of local fully qualified paths of Firehose GISTIC2 CNA file(s).
#'@param dCNA The local fully qualified path of a directory to store the rCNA RDS file(s).
#'@param silent Boolean. Default: FALSE. Whether or not to write progress information to console.
#'@param noLog Boolean. Default: FALSE. Whether or not to log results.
#
# Opens all files in a requested directory,
# reads them into a data frame, and save the resulting object in RDS format.
#
#'@family importM.COSMIC, importM.TCGA, importM.GISTIC2, importM.ProjectGenie
#
#'@example
# fFHG <- "../inst/extdata/devCNA.txt"
# dCNA <- "../inst/extdata/dCNA"
# silent <- FALSE
# noLog <- FALSE
#
# importCNA.GISTIC2(fFHG, dCNA, silent = FALSE, noLog = FALSE)
#
# Input - txt (table)
# first 3 columns are gene, NIH locus ID and cytoband, remaining columns are sample numbers
# each row for one gene
# copy number values are in table
#
# Output
# Rownames: HGNC gene symbols.
# Column names: Tumor Sample Barcodes
# Values: CNA as copynumber -2. Floating point number or NA
#
# TODO
# MAKE SURE PATH GIVEN IS VALID AND FILE IS NOT EMPTY/CORRUPT
# READ
# GET RID OF INFORMATION NOT NEEDED (FILTER)
# SAVE AS DATAFRAME AS RDS,  AND STORE IN DIR
# PRINT STATUS
