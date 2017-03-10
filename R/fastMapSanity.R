# fastMapSanity.R

#' Check if given input is a valid key or valid value
#'
#' \code{fastMapSanity} Will check if \code{toCheck} contains an illegal
#' character (Allowed characters: alphanumeric, -, _).
#'
#' @param toCheck Given input to check.
#'
#' @family fastMap functions
#'
#' @examples
#' fastMapSanity("ENSP00000245105")
#' fastMapSanity("A2M")
#' @export
fastMapSanity <- function(toCheck) {
    if (!grepl("^[a-zA-Z0-9_-]+$", toCheck)) {
        return(FALSE)
    }
    return(TRUE)
}

# [END]
