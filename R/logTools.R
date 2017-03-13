# logTools.R
#
# Utility functions logging and object metadata
#

#' Define the log file name
#'
#' \code{logFileName} returns path and name of a file for logging events, or
#'   sets the global option rete.logfile.
#'
#'   If fPath is missing, getwd() will be used as the path. If fPath ends with
#'   a "/", that character is removed since the file.path() function will add
#'   a "/" which would result in "//" if one already exists. If "//" is the
#'   intended path, add an extra "/". The default file
#'   name is structured as "rete_YYY-MM-DD.d.log" where YYY-MM-DD is the value
#'   of Sys.Date() and d is an integer that
#'   is one larger than the highest integer in a filename structured this way
#'   in the fPath directory. Thus the file name is unique. If the requested
#'   filename does not end with the extension ".log", that extension is added.
#'   If you must have a logfile named with a different extension, you can set
#'   it directly via options("rete.logfile") <- <my.special.name>.
#'   If setOption is TRUE, the function sets the global option rete.logfile.
#'
#' @param fPath A path to a log-file directory. If missing the path is set to
#'   getwd().
#' @param fName Name of file to which log messages will be appended. If
#'   missing, a new filename is generated.
#' @param setOption Whether to set the global option rete.logfile.
#'   Defaults to FALSE.
#' @return Path and name of a file that can be used to log events.
#'
#' ## @family TBD
#'
#'   ## @seealso \code{\link{logMessage}} ...
#'
#' @examples
#' logFileName()
#' \dontrun{
#' logFileName(fPath = "./logs", setOption = TRUE)
#' }
#' @export
logFileName <- function(fPath = getwd(), fName, setOption = FALSE) {

    # === VALIDATE fPath =======================================================
    r <- .checkArgs(fPath, like = "DIR", checkSize = TRUE)
    if(length(r) > 0) {
        stop(r)
    }

    # ensure path does not end with a "/" because we later use file.path()
    # and that adds a "/"
    fPath <- gsub("/$", "", fPath)

    if (missing(fName) || fName == "") {
        today <- Sys.Date()
        files <- list.files(path = fPath,
                            pattern = paste(today,
                                            "\\.[0-9]+\\.log$",
                                            sep = ""),
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

    # === VALIDATE fName, setOption ============================================
    r <- c(r, .checkArgs(fName, like = "a", checkSize = TRUE))
    r <- c(r, .checkArgs(setOption, like = TRUE, checkSize = TRUE))
    if(length(r) > 0) {
        stop(r)
    }

    # === add ".log" extension if there is none ================================

    if (! grepl("\\.log$", fName)) {
        fName <- paste(fName, ".log", sep = "")
    }

    if(setOption) {
        options(rete.logfile = file.path(fPath, fName))
    }
    return(file.path(fPath, fName))
}


#' Write a message to the current log file.
#'
#' \code{logMessage} uses cat() to append a message to the current log
#'   file. The file name is taken from the rete.logfile global option.
#'
#' The function will stop() if message is not of mode, type and class character.
#' If a line is not terminated with a linebreak, a linebreak is appended.
#'
#' @param message a character object or vector of character objects.
#' @return N/A. message is appended to the current logfile.
#'
#' @family log file functions
#'
#'   ## @seealso \code{\link{logFileName}}
#'
#' @examples
#' \dontrun{
#'   mesg <- c("note > ", Sys.Date(), "Re-run analysis.")
#'   logMessage(mesg)
#' }
#' @export
logMessage <- function(message) {

    checkReport <- .checkArgs(message, like = character())
    if(length(checkReport) > 0) {
        stop(checkReport)
    }

    # append "\n" to lines that don't have it.
    for (i in grep("\n$", message, invert = TRUE)) {
        message[i] <- paste(message[i], "\n", sep = "")
    }

    cat(message,
        file = unlist(getOption("rete.logfile")),
        sep = "",
        append = TRUE)

}

#' Attach a UUID to an object
#'
#' \code{attachUUID} uses the uuid package to attach a $UUID attribute to
#' a given object.
#'
#' The overwrite flag can be understood as follows:
#' - if the overwrite flag is TRUE, then whether or not the object had
#' a UUID previously or not, a new UUID attribute is attached
#' - if the overwrite flag is FALSE, then no new UUID attribute is attached
#'
#' @param object any R object
#' @param overwrite flag for indicating whether a UUID should be overwritten or not
#' @return NULL if no new UUID is attached, and the value of the attached UUID if a
#' new UUID was attached
#'
#' @family log file functions
#'
#'   ## @seealso \code{\link{logFileName}}
#'   ## @seealso \code{\link{logMessage}}
#'   ## @seealso \code{\link{extractAttributes}}
#'   ## @seealso \code{\link{logEvent}}
#'   ## @seealso \code{\link{findUUID}}
#'   ## @seealso \code{\link{getProvenance}}
#'
#' @examples
#' \dontrun{
#'   object <- "c"
#'   attachUUID(object)
#' }
#' @export
attachUUID <- function(object, overwrite = TRUE) {

    if (is.null(object)) {
        stop(object)
    }

    if (overwrite || is.null(attr(object, "UUID"))) {
        generatedUUID <- uuid::UUIDgenerate()
        attr(object, "UUID") <- generatedUUID
        return(object)
    } else {
        return(object)
    }

}

#' Extract attributes from object
#'
#' \code{extractAttributes} returns all of the attributes that are attached
#' to an object, formatted for logging
#'
#' Format of the extracted attributes:
#' attribute | <attrName1> | <attrValue1>
#' attribute | <attrName2> | <attrValue2>
#' attribute | <attrName3> | <attrValue3>
#'
#' @param object any R object
#' @return NULL if no attributes are extracted, text of attributes if attributes exist
#'
#' @family log file functions
#'
#'   ## @seealso \code{\link{logFileName}}
#'   ## @seealso \code{\link{logMessage}}
#'   ## @seealso \code{\link{attachUUID}}
#'   ## @seealso \code{\link{logEvent}}
#'   ## @seealso \code{\link{findUUID}}
#'   ## @seealso \code{\link{getProvenance}}
#'
#' @examples
#' \dontrun{
#'   object <- "c"
#'   attr(object, "UUID") <- uuid::UUIDgenerate()
#'   extractAttributes(object)
#' }
#' @export
extractAttributes <- function(object) {

    if (is.null(object)) {
        stop(object)
    }

    result <- ""
    for (name in names(attributes(object))) {
        attribute <- sub("\\s+$", "", sub("^\\s+", "", name))
        value <- sub("\\s+$", "", sub("^\\s+", "", attr(object, name)))

        result <- paste(result, "attribute\t|\t", name, "\t|\t", attr(object, name), "\n", sep="")
    }

    if (result == "") {
        return(NULL)
    } else {
        return(result)
    }

}

# [END]
