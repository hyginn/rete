# importSNV.TCGA.R

#' Import MAF files from source and converts them to rSNV files.
#'
#' \code{anuvar} Adds annotation for mutations in a SNV table from GeneData table.
#'
#' @section Validations:
#'
#' @section Add Annotation:
#'
#' @section Write Log: log section writes file processing progress along with other errors
#' that may occur.
#'
#' @param SNV table of rows containing mutations per sample.
#' @param geneData table containing a list of gene annotations.
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#'
#' @family
#'
#'   ## @seealso \code{\link{anuvar}} anuvar()
#'
#'   ## @examples ## \dontrun {
#'   ## anuvar(SNV, geneData) ## }
#' @export
#' @return SNV table with annotations added from the GeneData table.

anuvar <- function(          SNV,
                             geneData,
                             silent = FALSE,
                             writeLog = TRUE) {


    myNotes <- character()
    myCall <- character()

    # ==== PARAMETERS ==========================================================

    NL <- .PlatformLineBreak()
    # # Read header and one contents line
    # tmp <- readLines(fName, n = 2)

    # ==== VALIDATIONS =========================================================

    # General parameter checks
    cR <- character()
    cR <- c(cR, .checkArgs(SNV,          like = "a",        checkSize = TRUE))
    cR <- c(cR, .checkArgs(geneData,     like = "a",        checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,     like = logical(1), checkSize = TRUE))

    if(length(cR) > 0) {
        stop(cR)
    }

    # ==== Validata GeneData and SNV  =======================================

    # Validate that SNV table contatins valid column names
    # Validate GeneData table contains valid column names

    # ==== READ DATA ===========================================================

    # 1. FIND rows in SNV mutation table that exists in GeneData table

    # 2. ADD extra information from the GeneData table to SNV table
    # # exonStart - chromosomal coordinates of all exon first nucleotides for the canonical cDNA sequence
    # # exonEnd - chromosomal coordinates of all exon last nucleotides for the canonical cDNA sequence
    # # map - function to map a chromosome coordinate to a cDNA position for this gene
    # # cDNA - canonical cDNA sequence of the gene
    # # phastCons - phastCons conservation scores (* 1000) for each cDNA nucleotide (from UCSC)
    # # phyloP - phyloP conservation scores (* 1000) for each cDNA nucleotide (from UCSC)
    # # 3D map - rownames are the PDB-ID/chain, columns are the first and last residue for which 3D coordinates have been mapped.
    # # 3D - list indices are PDB-ID/chain. For each PDB-ID/chain there is a data frame in which rownames a
    # #      are the gene's sequence positions for which 3D coordinates have been mapped to this structure.
    # #      Columns are the the structure residue ID, amino acid type, and the x, y, and z coordinates of
    # #      the Cα atom and of the side-chain centroid.

    # 3. Validate Class using function ValidateMapping

    # ==== WRITE LOG ===========================================================

    if(writeLog) {

        myTitle <- "anuvar"

        # Compile function call record{
        myCall[1] <- "anuvar("
        myCall[2] <- sprintf("SNV = %s, ", SNV)
        myCall[3] <- sprintf("geneData = %s, ", geneData)
        myCall[4] <- sprintf("silent = %s, ", as.character(silent))
        myCall[5] <- sprintf("writeLog = %s, ", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")

        # Record progress information
        myNotes <- c(myNotes, paste("Number of annotations added -", length(fMAF)))

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
                 notes = myNotes)
    }
}

# [END]
