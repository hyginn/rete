# logTools.R
#
# Utility functions logging and object metadata
#


# ==== logFileName() ===========================================================

#' Define the log file name
#'
#' \code{logFileName} returns path and name of a file for logging events, or
#'   sets the global option rete.logfile.
#'
#'   If fPath is missing, the function checks whether getwd() contains a
#'   directory called "logs". If yes, log files will be written to that
#'   directory, if no, getwd() will be used as the path. This default behaviour
#'   also defines the initial log file path and name on loading the package. A
#'   single "/" (or "\"), path separator is removed from fPath if present, since
#'   the file.path() function will add a path separator which would result in
#'   doubling it otherwise. The default file name is structured as
#'   "rete_YYY-MM-DD.d.log" where YYY-MM-DD is the value of Sys.Date() and d is
#'   an integer that is one larger than the highest integer in a filename
#'   structured this way in the fPath directory. Thus the file name is unique.
#'   If the requested filename does not end with the extension ".log", that
#'   extension is added. If you must have a logfile named with a different
#'   extension, you can set it directly via options("rete.logfile") <-
#'   <my.special.name>. If setOption is TRUE, the function sets the global
#'   option rete.logfile to which \code{\link{logMessage}} appends log events.
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
logFileName <- function(fPath, fName, setOption = FALSE) {

    # === SET fPath DEFAULT ====================================================
    if (missing(fPath)) {
        if (dir.exists(file.path(getwd(), "logs"))) {
            # "logs" subdirectory exists ...
            fPath <- file.path(getwd(), "logs")
        } else {
            # "logs" subdirectory does not exist ...
            fPath <- getwd()
        }
    }

    # === VALIDATE fPath =======================================================
    r <- .checkArgs(fPath, like = "DIR", checkSize = TRUE)
    if(length(r) > 0) {
        stop(r)
    }

    # ensure path does not end with a "/" or "\" because we later use
    # file.path() and that adds a separator
    fPath <- gsub("[/\\]$", "", fPath)

    # === SET UNIQUE fName DEFAULT =============================================
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

    # === SET GLOBAL OPTION IF REQUESTED =======================================
    if(setOption) {
        options(rete.logfile = file.path(fPath, fName))
    }


    return(file.path(fPath, fName))
}


# ==== logMessage() ============================================================


#' Write a message to the current log file.
#'
#' \code{logMessage} uses cat() to append a message to the current log
#'   file. The file name is taken from the rete.logfile global option.
#'
#' The function will stop() if message is not of mode, type and class character.
#' On windows systems, cat() replaces \\n linebreaks with \\r\\n. Therefore
#' logMessage() converts all linebreaks internally to \\n before handing them to
#' cat().
#'
#' @param msg a character object or vector of character objects.
#' @return N/A. This function is invoked for its side-effect to appended text to
#'   the current logfile.
#'
#' @family log file functions
#'
#'   ## @seealso \code{\link{logFileName}}
#'
#' @examples
#' \dontrun{
#'   msg <- c("note > ", Sys.Date(), "Re-run analysis.")
#'   logMessage(msg)
#' }
#' @export
logMessage <- function(msg) {

    # Check that argument is of mode character
    r <- .checkArgs(msg, like = character())
    if(length(r) > 0) { stop(r) }

    # Remove all line-end linebreaks
    msg <- gsub("[\r\n]+$", "", msg)

    # Remove all \r characters to convert windows into unix linebreaks
    msg <- gsub("\r", "", msg)

    # Add \n to line-ends
    msg <- gsub("$", "\n", msg)

    # Append msg to log file
    cat(msg,
        file = unlist(getOption("rete.logfile")),
        sep = "",
        append = TRUE)

}


# ==== getUUID() ===========================================================

#' get a UUID to attach to an object
#'
#' \code{getUUID} uses the uuid package to prepare a UUID suitable to
#'                be attached to an object identified by objectName.
#'
#' The overwrite flag can be understood as follows:
#' - if the object does not have a $UUID attribute, a UUID is returned
#' - if the object has a $UUID attribute and "overwrite" is TRUE
#'   a new UUID is returned.
#' - if the object has a $UUID attribute and "overwrite" is FALSE, the
#'   old UUID is returned.
#'
#'   Note: no check is done whether the original $UUID attribute
#'   is a valid UUID.
#'
#' @param objectName char. Name of an R object
#' @param overwrite logical. Flag to control whether an existing UUID should be
#'                  overwritten. Default FALSE.
#' @return the value of the original UUID attribute, or a new UUID
#'
#' @family log file functions
#'
#'   ## @seealso \code{\link{logFileName}}
#'   ## @seealso \code{\link{logMessage}}
#'   ## @seealso \code{\link{logEvent}}
#'   ## @seealso \code{\link{findUUID}}
#'   ## @seealso \code{\link{getProvenance}}
#'
#' @examples
#' \dontrun{
#'   tmp <- "c"
#'   getUUID("tmp")
#' }
#' @export
getUUID <- function(objectName, overwrite = FALSE) {

    found <- FALSE
    for (frame in rev(sys.parents())) {
        if (exists(objectName, frame = frame, inherits = FALSE)) {
            obj <- get(objectName, envir = sys.frame(frame))
            found <- TRUE
            break()
        }
    }
    if (! found) { stop(sprintf("Object \"%s\" does not exist.", objectName)) }

    if (is.null(obj)) {
        stop(sprintf("Object \"%s\" is NULL.",
                     objectName))
    }


    # if (! is.null(get0(objectName, envir = parent.frame(99)))) {
    #     obj <- get(objectName, envir = parent.frame(99))
    # } else {
    #     obj <- get(objectName)
    # }


    if (is.null(attr(obj, "UUID")) || overwrite) {
        UUID <- uuid::UUIDgenerate()
    } else {
        UUID <- attr(obj, "UUID")
    }
    return(UUID)

}

# ==== .extractAttributes() ====================================================

.extractAttributes <- function(obName, role) {

    # Non-exported helper function returns all of the attributes that are
    # attached to an object identified by objectName, formatted for logging.
    # Function looks for the object first in the parent.frame(), then on the
    # normal search path. This fixes a problem that object could not be found if
    # the function was called via the test_that environment.

    # Format of the extracted attributes:
    # event | input  | attribute | <attrNameN> | <attrValueN>
    # ... or
    # event | output | attribute | <attrNameN> | <attrValueN>
    #
    # Parameters:
    #    obName: char. Name of a non-NULL R object
    #    role:   <input|output|using> signifies the object role in the calling
    #               function's workflow.
    # Value:
    #    A zero-length string if no attributes are extracted,
    #    a vector of formatted attribute descriptions otherwise.

    supportedRoles <- c("input", "output", "using")

    # ==== PARAMETER CHECKS ====================================================
    if (missing(role)) {
        stop("role parameter must be provided.")
    }

    r <- character()
    r <- c(r, .checkArgs(obName, like = "a", checkSize = TRUE))
    r <- c(r, .checkArgs(role,   like = "a", checkSize = TRUE))
    if(length(r) > 0) {
        stop(r)
    }

    if (! role %in% supportedRoles) {
        stop(sprintf("Expecting role from (\"%s\"), but got \"%s\".",
                     paste(supportedRoles, collapse = "\", \""), role))
    }
    # ==== GET OBJECT ==========================================================
    found <- FALSE
    for (frame in rev(sys.parents())) {
        if (exists(obName, frame = frame, inherits = FALSE)) {
            myOb <- get(obName, envir = sys.frame(frame))
            found <- TRUE
            break()
        }
    }
    if (! found) { stop(sprintf("Object \"%s\" does not exist.", obName)) }

    if (is.null(myOb)) {
        stop(sprintf("Object \"%s\" is NULL.", obName))
    }


    # ==== COMPILE OUTPUT STRING ===============================================
    result <- character()
    for (name in names(attributes(myOb))) {

        if (length(attr(myOb, name)) == 1) {
            att <- paste("\"",
                         attr(myOb, name),
                         "\"",
                         sep = "")
        } else if (length(attr(myOb, name)) <= 3) {
            att <- paste("(",
                         paste(attr(myOb, name), collapse = ", "),
                         ")",
                         sep = "")
        } else {
            att <- paste("(",
                         paste(attr(myOb, name)[1:3], collapse = ", "),
                         ", ... (",
                         length(attr(myOb, name)),
                         ") )",
                         sep = "")
        }

        result <- c(result,
                    sprintf("event | %-6s | attribute | %-12s | %s",
                            role, name, att))
    }

    return(result)

}

# ==== logEvent() ==============================================================

#' Formats an event description to be sent for attaching to the log file.
#'
#' \code{logEvent} generates an event message with:
#' - the event title
#' - the event call
#' - attributes of the input object
#' - attributes of the output object
#' - date and time of the event
#' - function version (TODO)
#' - file hashes of input (TODO)
#'
#' Format of the event description is as follows:
#' event | title  | <eventTitle>
#' event | time   | <eventDateTime>
#' event | call   | <eventCall>
#' (ToDo: event | functionVersion | <functionVersion>)
#' event | input  | attribute | <inputAttrName1> | <inputAttrValue1>
#' event | output | attribute | <outputAttrName1> | <outputAttrValue1>
#' event | end
#'
#'
#' The event message is terminated by an "end event" marker, and a blank
#' line, and it is handed off to logMessage() to append it to the
#' log file referenced in the global variable options("rete.logfile")
#'
#'
#' @param eventTitle the title of the event, from a
#'                   predetermined categorization of events
#' @param eventCall the function call in which the event occurred
#' @param input vector of object names with an "input" role in the calling
#'              function's workflow
#' @param notes a vector of strings comprising text to be
#'              incorporated into the event description.
#' @param output vector of object names  with an "output" role in the calling
#'               function's workflow
#' @return N/A A message is appended to the log file
#'             via \code{\link{logMessage}}.
#'
#' @family log file functions
#'
#'   ## @seealso \code{\link{logFileName}}
#'   ## @seealso \code{\link{logMessage}}
#'   ## @seealso \code{\link{getUUID}}
#'   ## @seealso \code{\link{findUUID}}
#'   ## @seealso \code{\link{getProvenance}}
#'
#' @examples
#' \dontrun{
#'     logEvent(eventTitle = "Test event",
#'              eventCall = "function1(test = \"test\")")
#' }
#' @export
logEvent <- function(eventTitle,
                     eventCall,
                     input = character(),
                     notes = character(),
                     output = character()) {

    # check parameters
    if (missing(eventTitle)) {
        stop("eventTitle is missing with no default.")
    }

    if (missing(eventCall)) {
        stop("eventCall is missing with no default.")
    }

    # validate parameter structure
    r <- character()
    r <- c(r, .checkArgs(eventTitle, like = "a", checkSize = TRUE))
    r <- c(r, .checkArgs(eventCall,  like = "a", checkSize = TRUE))
    r <- c(r, .checkArgs(input, like = character()))
    if (length(notes) > 0) {
        r <- c(r, .checkArgs(notes, like = character()))
    }
    r <- c(r, .checkArgs(output, like = character()))
    if(length(r) > 0) {
            stop(r)
    }

    # validate objects
    for (obName in c(input, output)) {
        # find the object
        for (frame in rev(sys.parents())) {
            if (exists(obName, frame = frame, inherits = FALSE)) {
                myOb <- get(obName, envir = sys.frame(frame))
                found <- TRUE
                break()
            }
        }
        if (! found) { stop(sprintf("Object \"%s\" does not exist.", obName)) }

        # validate the object
        if (is.null(myOb)) {
            stop(sprintf("Object \"%s\" is NULL.", obName))
        }
    }

    msg <- character()
    msg <- c(msg, paste("event | title  | ", eventTitle, sep = ""))
    msg <- c(msg, paste("event | time   | ", Sys.time(), sep = ""))
    msg <- c(msg, paste("event | call   | ", eventCall, sep = ""))

    # main event logging
    for (i in input) {
        msg <- c(msg, sprintf("event | input  | \"%s\"", i))
        msg <- c(msg, .extractAttributes(i, "input"))
    }

    for (note in notes) {
        msg <- c(msg, sprintf("event | note   | %s", note))
    }

    for (o in output) {
        msg <- c(msg, sprintf("event | output | \"%s\"", o))
        msg <- c(msg, .extractAttributes(o, "output"))
    }

    msg <- c(msg, "event | end") # attach end marker
    msg <- c(msg, "") # attach an empty line for better readability

    logMessage(msg)
}


# ==== findUUID() ==============================================================

#' Finds a specific UUID in logged events
#'
#' \code{findUUID} returns the event blocks that contains mentions of a
#' requested UUID from the log files found in a given directory.
#'
#' @param uuid a universally unique identifier (UUID). Checked for valid format.
#' @param logDir the directory path in which to search the logfiles.
#'         Defaults to the path component of getOption("rete.logfile").
#' @param ext file extension regex pattern. Defaults to "\\.log$".
#' @param recursive Whether to descend into subdirectories. Defaults
#'                  to TRUE.
#' @return character() if uuid is not found in the log files.
#'         Otherwise, a vector with one complete event block per element, which
#'         has the filename that conteined the event attaches as a names()
#'         attribute.
#'
#' @family log file functions
#'
#'   ## @seealso \code{\link{logFileName}}
#'   ## @seealso \code{\link{logMessage}}
#'   ## @seealso \code{\link{getUUID}}
#'   ## @seealso \code{\link{logEvent}}
#'   ## @seealso \code{\link{getProvenance}}
#'
#' @examples
#' \dontrun{
#' NULL
#' }
#' @export
findUUID <- function(uuid,
                     logDir,
                     ext = "\\.log$",
                     recursive = TRUE) {

    if (missing(logDir)) {
        logDir <- dirname(unlist(getOption("rete.logfile")))
    }
    # check parameters
    r <- character()
    r <- c(r, .checkArgs(uuid, like = "UUID", checkSize = TRUE))
    r <- c(r, .checkArgs(logDir, like = "DIR", checkSize = TRUE))
    r <- c(r, .checkArgs(ext, like = "a", checkSize = TRUE))
    r <- c(r, .checkArgs(recursive, like = TRUE, checkSize = TRUE))
    if (length(r) > 0) {
        stop(r)
    }

    allLogs <- list.files(path = logDir,
                          pattern = ext,
                          full.names = TRUE,
                          recursive = recursive)

    if (length(allLogs) == 0) {
        stop(sprintf("No logfiles found in %s.", logDir))
    }

    allEvents <- character()

    for (fn in allLogs) {
        # read entire file into one string, then strsplit on the
        # boundaries of an event block
        tmp <- unlist(strsplit(readChar(fn, file.info(fn)$size),
                    "(?=\\nevent \\| title)|(?<=event \\| end)",
                    perl = TRUE))
        # find events that contain the UUID
        tmp <- tmp[grepl(uuid, tmp)]
        # attach the filename as name attribute
        names(tmp) <- rep(fn, length(tmp))
        allEvents <- c(allEvents, tmp)
    }

    # add two platform-appropriate linebreaks to each event block
    if (length(allEvents) > 0) {
        NL2 <- paste0(.PlatformLineBreak(), .PlatformLineBreak())
        for (i in 1:length(allEvents)) {
            allEvents[i] <- gsub("$", NL2, allEvents[i])
        }
    }

    return(allEvents)
}


#' Reconstructs provenance of object
#'
#' \code{getProvenance} reconstructs the provenance of an object back via its
#'                      intermediaries,  to all input files that contributed
#'                      to it.
#'
#' Format of the output... TBD.
#'
#'
#' @param uuid a universally unique identifier (UUID). Checked for valid format.
#' @param logDir the directory path in which to search the logfiles.
#'         Defaults to the path component of getOption("rete.logfile").
#' @param ext file extension regex pattern. Defaults to "\\.log$".
#' @param recursive Whether to descend into subdirectories. Defaults
#'                  to TRUE.
#' @return Informative message if uuid is not found in the log files.
#'         Otherwise, a structured report that describes the provenance of
#'         the object that is identified by uuid.
#'
#' @family log file functions
#'
#'   ## @seealso \code{\link{logFileName}}
#'   ## @seealso \code{\link{logMessage}}
#'   ## @seealso \code{\link{attachUUID}}
#'   ## @seealso \code{\link{findUUID}}
#'   ## @seealso \code{\link{logEvent}}
#'
#' @examples
#' \dontrun{
#' NULL
#' }
#' @export
getProvenance <- function(uuid,
                          logDir,
                          ext = "\\.log$",
                          recursive = TRUE) {



        # if (missing(logDir)) {
        #     logDir <- dirname(unlist(getOption("rete.logfile")))
        # }
        # # check parameters
        # r <- character()
        # r <- c(r, .checkArgs(uuid, like = "UUID", checkSize = TRUE))
        # r <- c(r, .checkArgs(logDir, like = "DIR", checkSize = TRUE))
        # r <- c(r, .checkArgs(ext, like = "a", checkSize = TRUE))
        # r <- c(r, .checkArgs(recursive, like = TRUE, checkSize = TRUE))
        # if (length(r) > 0) {
        #     stop(r)
        # }
        #
        # mentions <- character()
        # # Try finding uuid. Report if it could not be found at all.
        # mentions <- c(mentions,
        #               findUUID(uuid = uuid,
        #                        logDir = logDir,
        #                        ext = ext,
        #                        recursive = recursive))
        #
        # if (length(mentions) == 0) {
        #     stop("Requested UUID was not found in the log files.")
        # }
        #
        # # Otherwise build a tree
        # #  - root is the event that has uuid as output. There must
        # #    only be one such node.
        # #  - children are events that are linked via (inputUUID -> outputUUID)
        # #  - leafs are events in which the outputUUID is created from file.
        # #      (Functions that create objects from file must write a
        # #       structured event message to that effect.)
        # #  Needs a tree list structure that identifies nodes, and children
        # #  Needs helper function to parse an event and insert it's information
        # #  into the tree list.
        #
        # # Traverse the tree depth first.
        # # summarize each node, (what happened), and its children (how
        # # many, what are they, what are their UUIDs)
        #
        # # Output report.

    return("Not yet implemented.")

}
# [END]
