# fastMapUpdate.R

#' Modify a fastMap hash table.
#'
#' \code{fastMapUpdate} Update, add, or delete entries in the fastMap
#' hash table. Adding and updating are the same; if a key does/doesn't
#' exist, will set the key's value to \code{value}. If a lookup is made
#' for a non-existant key, NULL will be returned; setting a key's value to
#' NULL is essentially deleting. Before inserting into the hash table, convert
#' all letters in \code{value} to upper case.
#'
#' TODO: Vectorize.
#'
#' @param hashTable The fastMap hash table.
#' @param key The key to update.
#' @param value The HGNC symbol to associate with \code{key}.
#'
#' @family fastMap functions
#'
#' @seealso \code{\link{fastMapGenerate}} on how the hash table is structured.
#' \code{\link{fastMapSanity}} on acceptable keys and values.
#'
#' @examples
#' \dontrun{
#' fastMapUpdate(fastMapUniprot, "P01023", "A2M")
#' fastMapUpdate(fastMapENSP, "ENSP00000245105", "a2m-as1")
#' }
#' @export
fastMapUpdate <- function(hashTable, key, value) {
    # Key validation
    if (is.character(key)) {
        keySanity <- fastMapSanity(key)
    } else {
        keySanity <- FALSE
    }
    if (!keySanity) {
        errorMessage <- "Invalid key."
        stop(errorMessage)
    }

    # Value validation
    if (is.character(value)) {
        valueSanity <- fastMapSanity(value)
    } else if (is.null(value)) {
        valueSanity <- TRUE
    } else {
        valueSanity <- FALSE
    }
    if (!valueSanity) {
        errorMessage <- "Invalid value."
        stop(errorMessage)
    }

    # Update the hash table
    if (is.null(value)) {
        hashTable[[key]] <- value
    } else {
        hashTable[[key]] <- toupper(value)
    }
}

# [END]
