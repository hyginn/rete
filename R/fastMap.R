# fastMap.R

#' Map an ID to its associated HUGO gene symbol.
#'
#' \code{fastMap} accesses the associated hash table
#' (\code{fastMapUniProt} for UniProt IDs or
#' \code{fastMapUniENSP} for Ensembl protein IDs) and return the
#' associated HUGO gene symbol. Unmapped IDs will be return as is and
#' will be stored in an global option \code{rete.unmapped}.
#'
#' @param ID An unmapped ID.
#' @param type The type of the unmapped ID.
#' @return The HUGO gene symbol of \code{ID}
#'
#' @family fastMap functions
#'
#' @examples
#' fastmap("ENSG00000121410", type = "ENSP")
#' fastmap(c("Q9NQ94", "P01023"))
#' @export
fastMap <- function(ID, type = "UniProt") {

    # Check if associated table exists. Terminate if it doesn't.
    tableName <- paste0("fastMap", type)
    if (!exists(tableName)) {
        errorMessage <- paste("Hash table:", tableName, "is not loaded")
        stop(errorMessage)
    }

    unmappedCount <- 0
    IDLength = length(ID)
    out <- character(length = IDLength)
    for (i in 1:IDLength) {
        if (type == "UniProt") {
            lookUp <- fastMapUniProt[[ID[i]]]
        } else {
            lookUp <- fastMapENSP[[ID[i]]]
        }

        if (is.null(lookUp)) {
            unmappedCount <- unmappedCount + 1
            # append to rete.unmapped option
            options(rete.unmapped = c(getOption("rete.unmapped"), ID[i]))
            out[i] <- ID[i]
        } else {
            out[i] <- lookUp
        }
    }

    if (unmappedCount) {
        warningMessage <- paste(unmappedCount, "of", IDLength,
                                "IDs could not be mapped, see getOptions('rete.unmapped') to list them all.")
        warning(warningMessage)
    }
    return(out)
}

# [END]
