# checkGeneSymbols.R

#' Check whether gene symbols given are valid HGNC gene symbols.
#'
#' \code{isGeneSymbol} Recieves a vector as input and checks whether the input
#' vector contains valid HGNC gene symbols (case insensitive) by returning a
#' vector of length input, containing TRUE for gene symbols and FALSE for others.
#'
#' @section Checking symbols
#'   The function calls fmatch on each symbol given and checks the return value
#'   of fmatch to see if it is a valid gene symbol. If fmatch returns a
#'   value/index, the symbol must be a valid gene symbol; if fmatch returns
#'   NA, the symbol must be invalid. The function uses is.na to check if the
#'   return value of fmatch is NA.
#'
#' @param input A vector that must be checked for valid gene symbols.
#' @return A vector of length input containing TRUE for valid HGNC
#' gene symbols and FALSE for others.
#'
#'
#' @seealso \code{\link{generateHGNCHash}} Generates a data frame containing
#' all valid HGNC gene symbols, calls fmatch on the data frame so it will
#' create a hash of the data frame to be used by checkGeneSymbols.R.
#'
#' @examples
#' isGeneSymbol(c("123", "234"))
#' isGeneSymbol(c("A1BG", "a1bg", "a1Bg", "A1bG))
#' isGeneSymbol(c("123", "A1CF", "234", "a1bg"))

library(microbenchmark)
a <- c("a1bg", "A1BG")
for (i in 1:1000){
    a <- append(a, as.character(i))
}

isGeneSymbol <- function(input) {

    # Check the validity of the input.
    if (is.null(input) || length(input) == 0 || typeof(input) != "character" ||
        mode(input) != "character" || class(input) != "character") {
            stop("Invalid input!")
    }

    # Check the input vector against the geneNames data frame.
    return(!is.na(fmatch(toupper(input), geneNames$Gene_Symbol)))
}

# [END]
