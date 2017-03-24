# generateHGNCHash.R

#' Generates a data frame containing all valid HGNC gene symbols, calls fmatch
#' on the data frame so it will create a hash of the data frame to be used by
#' checkGeneSymbols.R.
#'
#' @param fURL The URL to the file containing valid HGNC gene symbols.
#'
#' @seealso \code{\link{isGeneSymbol}} is used to check if a gene symbol is a
#' valid HGNC gene symbol that is contained in inst/extdata/HGNCSymbols.rds.

library(fastmatch)

# Open the connection/file.
## Change fURL to be the URL of the file containing HGNC gene symbols.
## The file is assumed to contain the symbol, status and locus group of the gene.
fURL <- "http://www.genenames.org/cgi-bin/download?col=gd_app_sym&col=gd_status&col=gd_locus_type&col=gd_locus_group&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&hgnc_dbtag=on&submit=submit"
HGNCFile <- url(fURL, open = "r")

# Read the header and vertorize.
header <- readLines(HGNCFile, n = 1)
headerVector <- unlist(strsplit(header, "\t"))

# Find the index of the columns that contain the gene symbol, status, and locus
# group of the genes.
symbolCol <- which(headerVector == "Approved Symbol")
statusCol <- which(headerVector == "Status")
locusCol <- which(headerVector == "Locus Group")

# vectors containing the (important) statuses and locus groups of interest.
impStatuses <- c("Approved")
impLocusGroups <- c("protein-coding gene", "other", "phenotype")

# Read line by line.
geneSymbols <- c()
while (length(entry <- readLines(HGNCFile, n = 1)) > 0) {
    entryVector <- strsplit(entry, "\t")

    # Check the status and locus group of the gene.
    status <- entryVector[[1]][statusCol]
    locusGroup <- entryVector[[1]][locusCol]
    if ((status %in% impStatuses) && (locusGroup %in% impLocusGroups)){
        # append the current gene symbol to the list of valid gene symbols.
        geneSymbols <- append(geneSymbols, entryVector[[1]][symbolCol])
    }

}

close(HGNCFile)

# Create a data frame with a single column containing the gene symbols.
geneNames <- data.frame(
    Gene_Symbol = geneSymbols,
    stringsAsFactors = FALSE)

# Call fmatch, so that it will create a hash of the data frame.
fmatch("1", geneNames)

# [END]
