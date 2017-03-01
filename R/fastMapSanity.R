# fastMapSanity.R

#' Check if given input is a valid key or valid value
#'
#' \code{fastMapSanity} If the \code{type} is key, this will check if
#' \code{toCheck} is a string and does not contain an illegal character
#' (allowed characters are alphanumeric, -, _). If the \code{type} is value,
#' this will check if \code{toCheck} is a string or NULL. If it is a string,
#' it will check if \code{toCheck} does not contain an illegal character.
#'
#' @param toCheck Given input to check.
#' @param type The type of \code{toCheck}
#'
#' @family fastMap functions
#'
#' @examples
#' fastMapSanity("ENSP00000245105", "key")
#' fastMapSanity("A2M", "value")
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
