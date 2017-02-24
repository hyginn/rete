# zzz.R
#
# Package startup and unload functions


# We need to define this function here, because .onLoad() uses it
.makeLogFileName <- function(fName) {
    # make a unique log-file name or
    # return the requested filename
    if (missing(fName)) {
        today <- Sys.Date()
        files <- list.files(pattern = paste(today, "\\.[0-9]+\\.log$", sep = ""),
                            all.files = TRUE)
        if (length(files) > 0) {
            # get highest version of today's existing files
            m <- regexec("\\.([0-9]+)\\.log$", files)
            num <- max(as.numeric(unlist(regmatches(files, m))[c(FALSE, TRUE)]))
        } else {
            num <- 0
        }
        fName <- paste("rete_", Sys.Date(), ".", num + 1, ".log", sep = "")
    }
    return(fName)
}


.onLoad <- function(libname, pkgname) {

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Failed to find required package \"igraph\".")
    }

    if (!requireNamespace("readr", quietly = TRUE)) {
        stop("Failed to find required package \"readr\".")
    }

    op <- options()

    op.rete <- list(
        rete.logfile = .makeLogFileName()
        )

    toset <- !(names(op.rete) %in% names(op))

    if(any(toset)) {
        options(op.rete[toset])
    }

    invisible()
}


.onAttach <- function(libname, pkgname) {
    m <- character()
    m[1] <- "\nWelcome to rete.\n"
    m[2] <- sprintf("  Workflow will be logged to \"%s\"\n",
                    getOption("rete.logfile"))
    m[3] <- "  or type ?logfile for instructions to set a different name.\n\n"

    packageStartupMessage(paste(m, collapse=""))
}


# .onUnload <- function(libname, pkgname) {
#
# }



# [END]
