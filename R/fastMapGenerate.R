# fastMapGenerate.R

#' Generate the fastMap hash tables.
#'
#' \code{fastMapGenerate} will parse through an HGNC file and create a hash
#' table. Keys are the IDs found in \code{unmappedColName} and associated values
#' are IDs found in \code{HGNCSymbolColName}. The type attribute of the hash
#' table will be set to \code{type}. Hash table can be saved and can be loaded
#' in future sessions. If there are conflicting keys or a tab delimited key, keep
#' the first instance of the key or only store the first key respectively.
#'
#' @param fName The path to the HGNC gene symbol data set.
#' @param HGNCSymbolColName The column name of the HGNC symbol.
#' @param unmappedColName The column name the unmapped symbol.
#' @param type The type of the gene database of unmappedColName.
#' @param saveHashTable Boolean whether to save the hash table as an rds.
#' @param outputName The path to save the RDS.
#'
#' @seealso \code{\link{fastMapSanity}} on acceptable keys and values.
#'
#' @family fastMap functions
#'
#' @examples
#' \dontrun{
#' fastMapGenerate("hgnc_complete_set.txt", "symbol", "uniprot_ids",
#'                 type = "UniProt", outputName = "fastMapUniProt.rds")
#' }
#' @export
fastMapGenerate <- function(fName, HGNCSymbolColName, unmappedColName,
                            type, saveHashTable = TRUE, outputName) {
    if (saveHashTable) {
        if (missing(outputName)) {
            errorMessage <- "saveHashTable is set to TRUE but outputName is not specified."
            stop(errorMessage)
        } else if (is.character(outputName)) {
            # Check if rds extension is in outputName, warn if it doesn't
            if (!grepl("\\.rds$", outputName)) {
                warningMessage <- "Supplied outputName does not contain a .rds file extension."
                warning(warningMessage)
            }
        } else {
            errorMessage <- "Supplied outputName can only be a string."
            stop(errorMessage)
        }
    }
    hashTable <- new.env(hash = TRUE)
    listedUnique <- new.env(hash = TRUE)

    HGNCFile <- file(fName, open = "r")

    # Read header
    tmp <- readLines(HGNCFile, n = 1)
    header <- unlist(strsplit(tmp, "\t"))

    # Locate column for HGNC symbol and unmapped column
    symbolCol <- which(header == HGNCSymbolColName)
    unmappedCol <- which(header == unmappedColName)

    # Check if columns have been found
    if (length(symbolCol) == 0) {
        errorMessage <- sprintf("Symbol column: %s, is not found.",
                                HGNCSymbolColName)
        stop(errorMessage)
    } else if (length(unmappedCol) == 0) {
        errorMessage <- sprintf("Unmapped column: %s, is not found.",
                                unmappedColName)
        stop(errorMessage)
    }

    totalCount <- 0
    nonUniqueKeysCount <- 0
    unsanitizedKeys <- 0
    unsanitizedValues <- 0
    # Read line by line
    while (length(entry <- readLines(HGNCFile, n = 1)) > 0) {
        entryVector <- strsplit(entry, "\t")
        key <- entryVector[[1]][unmappedCol]
        value <- entryVector[[1]][symbolCol]
        # Check if the key is a pipe delimited. Choose the first entry.
        # skip if there is no id for the specified type
        if (key == "") {
            next
        }

        # Check key sanity
        if (is.character(key)) {
            keySanity <- fastMapSanity(key)
            if (!keySanity) {
                # Check if it was pipe delimited
                key <- strsplit(value, "|", fixed = TRUE)[[1]][1]
                keySanity <- fastMapSanity(key)
                if (!keySanity) {
                    unsanitizedKeys <- unsanitizedKeys + 1
                }
            }
        } else {
            keySanity <- FALSE
        }

        # Check value sanity
        if (is.character(value)) {
            valueSanity <- fastMapSanity(value)
        } else if (is.null(value)) {
            valueSanity <- TRUE
        } else {
            valueSanity <- FALSE
        }

        if (!valueSanity) {
            unsanitizedValues <- unsanitizedValues + 1
        }
        if (!keySanity || !valueSanity) {
            next
        }
        # Handle non unique keys
        if (!is.null(hashTable[[key]])) {
            nonUniqueKeysCount <- nonUniqueKeysCount + 1
            # If we included the original non conflicted key
            if (is.null(listedUnique[[key]])) {
                nonUniqueKeysCount <- nonUniqueKeysCount + 1
                listedUnique[[key]] <- TRUE
                initialPair <- sprintf("%s->%s", key, hashTable[[key]])
                # insert first instance of the key with conflict
                options(rete.conflictingKeys =
                            c(getOption("rete.conflictingKeys"), initialPair))
            }
            # insert the other instance of they key with conflict
            pair <- sprintf("%s->%s", key, value)
            options(rete.conflictingKeys =
                        c(getOption("rete.conflictingKeys"), pair))
        } else {
            hashTable[[key]] <- value
        }
        totalCount <- totalCount + 1
    }

    close(HGNCFile)

    if (unsanitizedKeys) {
        warningMessage <- sprintf("%d of %d keys are unsanitized and skipped.",
                                  unsanitizedKeys, totalCount)
        warning(warningMessage)
    }

    if (unsanitizedValues) {
        warningMessage <- sprintf("%d of %d values are unsanitized and skipped.",
                                  unsanitizedValues, totalCount)
        warning(warningMessage)
    }

    if (nonUniqueKeysCount) {
        warningMessage <- sprintf("%d of %d keys are not unique. See getOption('rete.conflictingKeys') to list them all.",
                                  nonUniqueKeysCount, totalCount)
        warning(warningMessage)
    }

    # Insert metadata
    attr(hashTable, "type")  <- type
    if (saveHashTable) {
        saveRDS(hashTable, outputName)
    }
    return(hashTable)
}

# [END]
