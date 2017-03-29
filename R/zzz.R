# zzz.R
#
# Package startup and unload functions




.onLoad <- function(libname, pkgname) {

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Failed to find required package \"igraph\".")
    }

    if (!requireNamespace("readr", quietly = TRUE)) {
        stop("Failed to find required package \"readr\".")
    }

    if (!requireNamespace("uuid", quietly = TRUE)) {
        stop("Failed to find required package \"uuid\".")
    }

    # Make list of rete global parameters

    # filepath of logfile
    op.rete <- list(rete.logfile = logFileName() )

    # make a gG Prototype
    op.rete[["rete.gGprototype"]] <-
        importNet.STRING(system.file("extdata",
                                     "STRINGprototype.txt",
                                     package="rete"),
                         silent = TRUE,
                         writeLog = FALSE)

    toSet <- !(names(op.rete) %in% names(options()))

    if(any(toSet)) {
        options(op.rete[toSet])
    }

    invisible()
}


.onAttach <- function(libname, pkgname) {
    m <- character()
    m[1] <- "\nWelcome to rete.\n"
    m[2] <- sprintf("  Object details will be logged to\n    \"%s\"\n",
                    unlist(getOption("rete.logfile")))
    m[3] <- "  or type ?logfile for instructions to set a different name.\n\n"

    packageStartupMessage(paste(m, collapse=""))
}


# .onUnload <- function(libname, pkgname) {
#
# }



# [END]
