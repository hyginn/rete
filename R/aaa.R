# fastCheckFactory.R
#
# Utility function to create closures for ID validation

.fastCheckFactory <- function(symbols) {
    # Creates a closure to be used to validate whether elements of a
    #   character vector are present in the "symbols" data frame.
    # Parameters:
    #     symbols: a one column dataframe
    # Value:
    #     fCheck: a closure

    symbols$symbol[1] # force access

    fCheck <- function(x) {
        # Checks whether elements of x are present in the "symbols" data
        # frame which was loaded into this closure's environment.
        #
        # Returns a vector of logicals of the same length as x
        if (missing(x) || length(x) == 0) { return(logical()) }

        return(!is.na(fastmatch::fmatch(toupper(as.character(x)),
                                        symbols$symbol)))

    }

    return(fCheck)
}

# [END]
