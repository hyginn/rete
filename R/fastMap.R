# fastMap.R

#' Map a UniProt or Ensemble Protein ID to its associated HUGO gene symbol.
#'
#' \code{fastMap} return the associated HUGO gene symbol.
#'
#' @param ID an unmapped ID.
#' @param type the type of the unmapped ID.
#' @return The HUGO gene symbol of \code{ID}
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
            lookUp <- fastMapEnsp[[ID[i]]]
        }
        # TODO: Handle this case
        # To log or not to log
        if (is.null(lookUp)) {
            unmappedCount <- unmappedCount + 1
            out[i] <- "NA"
        } else {
            out[i] <- lookUp
        }
    }

    if (unmappedCount) {
        warningMessage <- paste(unmappedCount, "of", IDLength,
                                "IDs could not be mapped, see
                                getOptions(rete.unmapped) to list them all.")
        warning(warningMessage)
    }
    return(out)
}

# [END]
