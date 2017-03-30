#' isHGNCsymbol.R
#'
#' Check whether gene symbols given are valid HGNC gene symbols.
#'
#' \code{isHGNCsymbol} Checks whether the elements of the input vector are
#'     valid HGNC gene symbols (case insensitive) by comparing to a subset
#'     of existing gene symbols.
#'
#' The subset of symbols used here contains only approved gene
#'     symbols of the following locus types:
#'   * gene with protein product
#'   * immunoglobulin gene
#'   * protocadherin
#'   * T-cell receptor gene
#'   * RNA: long non-coding, micro, ribosomal, transfer, small nuclear
#'          and nucleolar, Y and vault
#'   * endogenous retrovirus
#'
#' This function is a closure that contains the HGNC symbol table in its environment. It is produced as part of the .onLoad() tasks. The supporting table is stored in extdata/HGNCsymbols.RDS. The script that was used to generate this table is in scripts/generateHGNCtable.R.
#'
#' Checking is done in a case-insensitive manner.
#'
#' @param x A character vector
#' @return A vector of logicals of length x that contains TRUE for every
#'           element that is present in the HGNC symbol table and FALSE
#'           for all others.
#'
#' @examples
#' isHGNCsymbol()                                   # logical()
#' isHGNCsymbol(NULL)                               # logical()
#' isHGNCsymbol(0)                                  # FALSE
#' isHGNCsymbol("A2M")                              # TRUE
#' isHGNCsymbol(c("123", "234"))                    # vectorized
#' isHGNCsymbol(c("A1BG", "a1bg", "a1Bg", "A1bG"))  # case insensitive
#' x <- c(NA, "A1CF", NULL, "a1bg")                 # length preserving:
#' length(x)                                        #    3
#' isHGNCsymbol(x)                                  #    FALSE, TRUE, TRUE
#'
#' @export

isHGNCsymbol <- function(x) {
    stop("This function must be overwritten by a closure factory in .onLoad()")
}

tmp <- readRDS(system.file("extdata",
                           "HGNCsymbols.RDS",
                           package="rete"))

isHGNCsymbol <- .fastCheckFactory(tmp)

rm(tmp)


# [END]
