# multiply.R

#' Multiply two numbers.
#'
#' \code{multiply} returns the product of its two arguments.
#'
#' This is an inconvenience function that achieves exactly the same as the
#' \code{"*"} operator, which it uses internally, but with more characters to
#' type.
#' @section QA warning: Actually using this function instead of "*" will get you
#'   in trouble with QA.
#'
#' @param a A number.
#' @param b A number.
#' @return The product of \code{x} and \code{y}. It will be cast to the most
#'   general type of its arguments.
#'
#' @family inconvenience functions
#'
#' @seealso \code{\link[base]{prod}} for products, \code{\link[base]{*}} for the
#'   arithmetic operator.
#'
#' @examples
#' multiply(1, 1)
#' multiply(1i, 1i)
#' multiply(13, 3)
multiply <- function(a, b) {
    return(a * b)
}

# [END]