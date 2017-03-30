# utils.R
#
# internal (non-exported) utility functions

.PlatformLineBreak <- function() {
    # returns "\r\n" on windows, "\n" otherwise
    if (.Platform$OS.type == "windows") {
        NL <- "\r\n"
    } else {
        NL <- "\n"
    }
    return(NL)
}


.pBar <- function(i, l, nCh = 50) {
    # Draw a progress bar in the console
    # i: the current iteration
    # l: the total number of iterations
    # nCh: width of the progress bar
    WID <- min(l-1, nCh)
    ticks <- round(seq(1, l-1, length.out = WID))
    if (i < l) {
        if (any(i == ticks)) {
            p <- which(i == ticks)
            p1 <- paste(rep("#", p), collapse = "")
            p2 <- paste(rep("-", WID - p), collapse = "")
            cat(sprintf("\r|%s%s|", p1, p2))
            utils::flush.console()
        }
    }
    else { # done
        cat("\n")
    }
}


.checkArgs <- function(x, like, checkSize = FALSE) {
    # Purpose:
    #     checks for argument mode, type, class and other attributes
    #
    # Parameters:
    #             x: the parameter that was given
    #          like: a prototype that x is being compared to
    #     checkSize: if TRUE
    #                -  check length for all 1D objects
    #                -  check dim for all kD, k>1 models
    #                -  check vcount(), ecount() for igraph graphs
    # Details:
    #     like can be a keyword that can trigger specific checks:
    #        DIR: string is an existing directory
    #        FILE_E: File exists
    #        FILE_W: File exists and is writable
    #        UUID: string is a valid UUID (not NIL)
    #
    # Value:
    #     report: "" if everything is oK
    #             a descriptive message if there is an error
    # Examples:
    #     report <- .checkArgs(ID, like = character())
    #     report <- .checkArgs(limits, like = c(1.0, 1.0), checkSize = TRUE)
    #     report <- .checkArgs(EGG, like = getOptions("rete.EGGprototype"))
    #
    # ToDo:
    #     - Check specifically for a versions graph attribute if this is
    #           an igraph object produced by a rete function.
    #     - Check individual elements of like if like is vector and
    #           checksize has passed.
    #     - Check size first.
    #     - Add parameter exhaustive = FALSE and internal flag ABORT to
    #           abort after first failed check.
    #

    name <- deparse(substitute(x)) # the name of the parameter to check

    report <- character()  # init
    NL <- .PlatformLineBreak()

    # === keyword specific checks ==============================================
    if (class(like) == "character" && length(like) == 1 &&
        class(x)    == "character" && length(x)    >  0) {
        for (el in x) {
            if (like == "DIR" && ! dir.exists(el)) {
                # el must be a valid directory
                report <- c(report,
                            sprintf(".checkArgs> \"%s\" %s%s%s%s",
                                    name,
                                    "error: directory \"",
                                    el,
                                    "\" does not exist.",
                                    NL))
            }
            if (like == "FILE_E" && ! file.exists(el)) {
                # el must be an existing file (or directory)
                report <- c(report,
                            sprintf(".checkArgs> \"%s\" %s%s%s%s",
                                    name,
                                    "error: file \"",
                                    el,
                                    "\" does not exist.",
                                    NL))
            }
            if (like == "FILE_W" && file.access(el, mode = 2) != 0) {
                # el must be a writable filename (i.e. must also exist)
                report <- c(report,
                            sprintf(".checkArgs> \"%s\" %s%s%s%s",
                                    name,
                                    "error: \"",
                                    el,
                                    "\" is not a writable file.",
                                    NL))
            }
            if (like == "UUID") {
                patt <- paste0("\\{?[0-9a-f]{8}-?", # braces, hyphens optional
                               "[0-9a-f]{4}-?",
                               "[1-5]",   # version
                               "[0-9a-f]{3}-?",
                               "[89ab]",  # variant
                               "[0-9a-f]{3}-?",
                               "[0-9a-f]{12}\\}?")
                # Note: this rejects the NIL UUID event though it is
                # syntactically correct as per RFC4122, in that it fails the
                # semantic requirements for this package.
                if (! grepl(patt, el, ignore.case = TRUE)) {
                    report <- c(report,
                                sprintf(".checkArgs> \"%s\" %s%s%s%s",
                                        name,
                                        "error: \"",
                                        el,
                                        "\" is not a valid UUID.",
                                        NL))
                }
            }
        }
    }

    # === mode/typeof/class checks =============================================
    if (mode(x) != mode(like)) {
        report <- c(report,
                    sprintf(".checkArgs> \"%s\" %s%s%s%s%s%s",
                            name,
                            "mode error: argument has mode \"",
                            mode(x),
                            "\" but function expects mode \"",
                            mode(like),
                            "\".",
                            NL))
    }

    if (typeof(x) != typeof(like)) {
        report <- c(report,
                    sprintf(".checkArgs> \"%s\" %s%s%s%s%s%s",
                            name,
                            "type error: argument has type \"",
                            typeof(x),
                            "\" but function expects type \"",
                            typeof(like),
                            "\".",
                            NL))
    }

    if (checkSize) {
        if (is.null(dim(like))) { # expect 1D objects
            if (!is.null(dim(x))) { # x is not 1D
                report <- c(report,
                            sprintf(".checkArgs> \"%s\" %s%s%s%s",
                                    name,
                                    "dimension error: argument has dim (",
                                    paste(dim(x), collapse = ", "),
                                    ") but function expects 1D object.",
                                    NL))

            } else if (length(x) != length(like)) { # x is 1D but unequal length
                report <- c(report,
                            sprintf(".checkArgs> \"%s\" %s%s%s%s%s%s",
                                    name,
                                    "length error: argument has length ",
                                    length(x),
                                    " but function expects length ",
                                    length(like),
                                    ".",
                                    NL))
            }
        } else { # expect xD objects
            if (is.null(dim(x)) || dim(x) != dim(like)) { # x is 1D or != dim()
                report <- c(report,
                            sprintf(".checkArgs> \"%s\" %s%s%s%s%s%s",
                                    name,
                                    "dimension error: argument has dim (",
                                    if (is.null(dim(x))) "NULL" else
                                        paste(dim(x), collapse = ", "),
                                    ") but function expects dim (",
                                    paste(dim(like), collapse = ", "),
                                    ").",
                                    NL))
            }
        }
    } # if (checkSize)

    if (class(x) != class(like)) {
        report <- c(report,
                    sprintf(".checkArgs> \"%s\" %s%s%s%s%s%s",
                            name,
                            "class error: argument has class \"",
                            class(x),
                            "\" but function expects class \"",
                            class(like),
                            "\".",
                            NL))
    }


    # === class igraph specific checks =========================================
    if (class(like) == "igraph") {
        if (length(igraph::graph_attr(x)) != length(igraph::graph_attr(like))) {
            report <- c(report,
                        sprintf(".checkArgs> \"%s\" %s%d%s%d%s%s",
                                name,
                                "graph attribute error: argument has ",
                                length(igraph::graph_attr(x)),
                                " attributes but function expects ",
                                length(igraph::graph_attr(like)),
                                ".",
                                NL))
        } else if (length(igraph::graph_attr(like)) > 0) { # attributes exist
            if (!isTRUE(all.equal(sort(names(igraph::graph_attr(x))),
                                  sort(names(igraph::graph_attr(like)))))) {
                report <- c(report,
                            sprintf(".checkArgs> \"%s\" %s%s%s%s%s%s%s",
                                    name,
                                    "graph attribute name error:",
                                    " argument has names (",
                                    paste(names(igraph::graph_attr(x)),
                                          collapse =", "),
                                    " but function expects (",
                                    paste(names(igraph::graph_attr(like)),
                                          collapse =", "),
                                    ").",
                                    NL))
            } else { # both objects' attribute names are equal
                # check graph object "version" attribute
                if ("version" %in% names(igraph::graph_attr(like)) &&
                    igraph::graph_attr(x)$version !=
                    igraph::graph_attr(like)$version) {
                    report <- c(report,
                                sprintf(".checkArgs> \"%s\" %s%s%s%s%s%s%s",
                                        name,
                                        "graph version error:",
                                        " argument has version \"",
                                        igraph::graph_attr(x)$version,
                                        "\" but function expects \"",
                                        igraph::graph_attr(like)$version,
                                        "\".",
                                        NL))
                }
            }
        }
    }

    return(report)
}

# [END]
