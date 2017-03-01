# fastMapUpdate.R

#' Modify a fastMap hash table.
#'
#' \code{fastMapUpdate} Update, add, or delete entries in the fastMap
#' hash table. Adding and updating are the same; if a key does/doesn't
#' exist, will set the key's value to \code{value}. If a lookup is made
#' for a non-existant key, NULL will be returned; setting a key's value to
#' NULL is essentially deleting.
#'
#' @param hashTable The fastMap hash table.
#' @param key The key to update.
#' @param value The value to insert to the associated key.
#'
#' @family fastMap functions
#'
#' @seealso \code{\link{fastMapGenerate}} on how the hash table is structured.
#' \code{\link{fastMapSanity}} on acceptable keys and values.
#'
#' @examples
#' fastMapUpdate(fastMapUniprot, "P01023", "A2M")
#' fastMapUpdate(fastMapENSP, "ENSP00000245105", "a2m-as1")
#' @export
fastMapUpdate <- function(hashTable, key, value) {
    # Key and value validation
    if (!fastMapSanity(key, "key")) {
        errorMessage <- "Invalid key."
        stop(errorMessage)
    } else if (!fastMapSanity(value, "value")) {
        errorMessage <- "Invalid value."
        stop(errorMessage)
    }
    if (is.null(value)) {
        hashTable[[key]] <- value
    } else {
        hashTable[[key]] <- toupper(value)
    }
}

# [END]
