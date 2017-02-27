# fastMapUpdate.R

#' Modify a fastMap hash table.
#'
#' \code{fastMapUpdate} Update, add, or delete entries in the fastMap
#' hash table. Adding and updating are the same; if a key does/doesn't
#' exist, will set the key's value to \code{value}. If a lookup is made
#' for a non-existant, NULL will be returned; setting a key's value to
#' NULL is essentially deleting. Sanity checks are included to only allow
#' acceptable characters for both \code{key} and \code{value}. \code{value}
#' can only be a string or NULL.
#'
#' @param hashTable The fastMap hash table.
#' @param key The key to update.
#' @param value The value to insert to the associated key.
#'
#' @family fastMap functions
#'
#' @seealso \code{\link{fastMapGenerate}} on how the hash table is structured.
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
        if (!grepl("^[a-zA-Z0-9_-]+$", key)) {
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
        if (!grepl("^[a-zA-Z0-9_-]+$", value)) {
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
