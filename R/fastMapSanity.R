# fastMapSanity.R

#' Check if given input is a valid key or valid value
#'
#' \code{fastMapUpdate} Checks if the key is a string and doesn not contain
#' an illegal character (allowed are alphanumeric, -, _). Check sif the value is
#' a string or NULL. If it is a string, it also checks if value does not contain
#' an illegl character.
#'
#' @param toCheck Given input to check.
#' @param type The type of \code{toCheck}
#'
#' @family fastMap functions
#'
#' @examples
#' fastMapSanity("ENSP00000245105", "key")
#' fastMapSAnity(""A2M", "value")
fastMapSanity <- function(toCheck, type) {
    # Check types of the input
    if (type == "key") {
        if (!is.character(toCheck)) {
            return(FALSE)
        }
    } else if (type == "value") {
        if (!is.character(toCheck) && !is.null(toCheck)) {
            return(FALSE)
        }
    } else {
        return(FALSE)
    }
    # Alphanumeric check for characters only
    if (is.character(toCheck)) {
        if (!grepl("^[a-zA-Z0-9_-]+$", toCheck)) {
            return(FALSE)
        }
    }
    return(TRUE)
}

# [END]
