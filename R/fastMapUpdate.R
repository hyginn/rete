# fastMapUpdate.R

#' Modify a fastMap hash table.
#'
#' \code{fastMapUpdate} Update, add, or delete entries in the fastMap hash table.
#'
#' @param hashTable the fastMap hash table.
#' @param key an unmapped ID.
#' @param value the mapped HUGO gene symbol to update.
#'
#' @family fastMap functions
#'
#' @examples
#' fastMapUpdate(fastMapUniprot, "P01023", "A2M")
#' fastMapUpdate(fastMapENSP, "ENSG00000245105", "a2m-as1")
#' @export
fastMapUpdate <- function(hashTable, key, value) {

    # Key validation
    if (is.character(key)) {
        # Check for invalid characters
        # Alphanumeric, _, - only
        if (!grepl("^[a-zA-Z0-9_-]+", key)) {
            errorMessage <- "Value contains an illegal character."
            stop(errorMessage)
        }
    } else {
        errorMessage <- "Keys can only be a string."
        stop(errorMessage)
    }

    # Value validation
    if (is.null(value)) {
        hashTable[[key]] <- value
    } else if (is.character(value)) {
        # Check for invalid characters
        # Alphanumeric, _, - only
        if (!grepl("^[a-zA-Z0-9_-]+", value)) {
            errorMessage <- "Value contains an illegal character."
            stop(errorMessage)
        }
        hashTable[[key]] <- toupper(value)
    } else {
        errorMessage <- "Value can only be a string or NULL."
        stop(errorMessage)
    }
}

# [END]
