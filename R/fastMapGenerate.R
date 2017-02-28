# fastMapGenerate.R

#' Generate the fastMap hash tables.
#'
#' \code{fastMapGenerate} will parse through an hgnc file and create a hash
#' table depending on \code{type}. Also saves the generated hash table as an
#' .rda which can be loaded in future sessions. If there conflicting keys, only
#' keep the first instance.
#'
#' @param fName The path to the HUGO gene symbol data set.
#' @param type The type to parse for.
#' @param saveHashTable Boolean whether to save the hash table as an rda.
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
    hashTable <- new.env(hash = TRUE)
    listedUnique <- new.env(hash = TRUE)

    hgncFile <- file(fName, open = "r")

    # Read header
    tmp <- readLines(hgncFile, n = 1)
    header <- unlist(strsplit(tmp, "\t"))

    # Locate col for symbol and associated ID
    symbolCol <- which(header == "symbol")
    unmappedCol <- which(header == unmappedColName)

    totalCount <- 0
    nonUniqueKeysCount <- 0
    # Read line by line
    while(length(entry <- readLines(hgncFile, n = 1)) > 0) {
        entryVector <- (strsplit(entry, "\t"))
        key <- entryVector[[1]][unmappedCol]
        # skip if there is no id for the specified type
        if (key == "") {
            next
        }
        value <- entryVector[[1]][symbolCol]
        # Handle non unique keys
        if (!is.null(hashTable[[key]])) {
            nonUniqueKeysCount <- nonUniqueKeysCount + 1
            # If we included the original non conflicted key
            if (is.null(listedUnique[[key]])) {
                nonUniqueKeysCount <- nonUniqueKeysCount + 1
                listedUnique[[key]] <- TRUE
                initialPair <- paste0(key, "->", hashTable[[key]])
                # insert first instance of the key with conflict
                options(rete.conflictingKeys =
                            c(getOption("rete.conflictingKeys"), initialPair))
            }
            # insert the other instance of they key with conflict
            pair <- paste0(key, "->", value)
            options(rete.conflictingKeys =
                        c(getOption("rete.conflictingKeys"), pair))
        } else {
            hashTable[[key]] <- value
        }
        totalCount <- totalCount + 1
    }

    close(hgncFile)

    if (nonUniqueKeysCount) {
        warningMessage <- paste(nonUniqueKeysCount, "of",
                                totalCount,
                                "keys are non-unique. See getOption('rete.conflictingKeys') to list them all")
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
