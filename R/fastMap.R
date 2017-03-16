# fastMap.R

#' Map an ID to its associated HGNC gene symbol.
#'
#' \code{fastMap} accesses the associated hash table
#' (\code{fastMapUniProt} for UniProt IDs or
#' \code{fastMapUniENSP} for Ensembl protein IDs) and returns the
#' associated HGNC gene symbol. Unmapped IDs will be returned as is and
#' will be stored in an global option \code{rete.unmapped}.
#'
#' @param ID An unmapped ID.
#' @param hashTable The hash table to perform lookups on.
#' @param type The type of the unmapped ID.
#' @return The mapped HGNC symbol or \code{ID} if not found.
#'
#' @family fastMap functions
#'
#' @examples
#' \dontrun{
#' fastMap("ENSG00000121410", fastMapENSP, type = "ENSP")
#' fastMap(c("Q9NQ94", "P01023"), fastMapUniProt)
#' }
#' @export
fastMap <- function(ID, hashTable, type = "UniProt") {
    if (!is.null(attributes(hashTable)$type)) {
        if (type != attributes(hashTable)$type) {
            errorMessage <- sprintf("Expected hash table type: %s. Supplied hash table type: %s",
                                    type, attributes(hashTable)$type)
            stop(errorMessage)
        }
    } else {
        errorMessage <- "Supplied hash table does not have a type attribute."
        stop(errorMessage)
    }
    IDLength = length(ID)
    out <- character(length = IDLength)
    for (i in 1:IDLength) {
        lookUp <- hashTable[[ID[i]]]
        if (is.null(lookUp)) {
            out[i] <- NA
        } else {
            out[i] <- lookUp
        }
    }

    # Address unmapped IDs
    unmappedCount <- sum(is.na(out))
    if (unmappedCount) {
        iNA <- which(is.na(out))
        out[iNA] <- ID[iNA]
        # Export unmapped IDs to the global option rete.unmapped
        options(rete.unmapped = ID[iNA])
        warningMessage <- sprintf("%d of %d IDs could not be mapped, see getOptions('rete.unmapped') to list them all.",
                                  unmappedCount, IDLength)
        warning(warningMessage)
    }
    return(out)
}

# [END]
