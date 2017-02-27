# fastMapGenerate.R

#' Generate the fastMap hash tables.
#'
#' \code{fastMapGenerate} will parse through an hgnc file
#' and create a hash table depending on \code{type}. Also
#' saves the generated hash table as an .rda which can be
#' loaded in future sessions.
#'
#' @param fName the path to the HUGO gene symbol data set.
#' @param type the type to parse for.
#' @param saveHT boolean whether to save the hash table as an rda.
#'
#' @family fastMap functions
#' @export
fastMapGenerate <- function(fName, type, saveHashTable=TRUE) {

    if (type == "UniProt") {
        unmappedColName <- "uniprot_ids"
    } else if (type == "ENSP") {
        unmappedColName <- "ensembl_gene_id"
    } else {
        errorMessage <- "Type can only be UniProt or ENSP"
        stop(errorMessage)
    }
    hashTable <- new.env(hash=TRUE)
    nonUnique <- new.env(hash=TRUE)

    hgncFile <- file(fName, open = "r")

    # Read header
    tmp <- readLines(hgncFile, n = 1)
    header <- unlist(strsplit(tmp, "\t"))

    # Locate col for symbol and associated ID
    symbolCol <- which(header == "symbol")
    unmappedCol <- which(header == unmappedColName)

    totalCount <- 0
    nonUniqueKeysCount <- 0
    nonUniqueValuesCount <- 0
    # Read line by line
    while(length(entry <- readLines(hgncFile, n = 1)) > 0) {
        entryVector <- (strsplit(entry, "\t"))
        key <- entryVector[[1]][unmappedCol]
        # skip if there is no id for the specified type
        if (key == "") {
            next
        }
        value <- entryVector[[1]][symbolCol]
        # TODO: Address non-unique keys
        if (!is.null(hashTable[[key]])) {
            nonUniqueKeysCount <- nonUniqueKeysCount + 1
        }
        # TODO: Address non-unique values
        if (!is.null(nonUnique[[value]])) {
            nonUniqueValuesCount <- nonUniqueValuesCount + 1
        }
        hashTable[[key]] <- value
        nonUnique[[value]] <- key
        totalCount <- totalCount + 1
    }

    if (nonUniqueKeysCount) {
        warningMessage <- paste(nonUniqueKeysCount, "keys out of", totalCount, "are the same.")
        warning(warningMessage)
    }

    if (nonUniqueValuesCount) {
        warningMessage <- paste(nonUniqueValuesCount, "values out of", totalCount, "are the same.")
        warning(warningMessage)
    }

    # Store hash table as a global variable and save an RDS.
    if (type == "UniProt") {
        fastMapUniProt <<- hashTable
        if (saveHashTable) {
            save(fastMapUniProt, file = "fastMapUniProt.rda")
        }
    } else {
        fastMapENSP <<- hashTable
        if (saveHashTable) {
            save(fastMapENSP, file = "fastMapENSP.rda")
        }
    }

}

# [END]
