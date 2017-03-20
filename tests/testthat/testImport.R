# testImport.R
#
#
context("import networks")

# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
# ==== END SETUP AND PREPARE ===================================================


# ==== importNet.STRING() ======================================================

# ToDo test handling of corrupt STRING files

test_that("importNet.STRING() rejects erroneous parameters", {
    fN <- "dataSTRING.txt"

    # ==== fName
    expect_error(
        gG <- importNet.STRING(fName = NULL,
                               net = "experimental",
                               dropUnmapped = FALSE,
                               silent = TRUE,
                               writeLog = FALSE),
        '.checkArgs> "fName" mode error')

    expect_error(
        gG <- importNet.STRING(fName = c("fee", "fie", "fo", "fum"),
                               net = "experimental",
                               dropUnmapped = FALSE,
                               silent = TRUE,
                               writeLog = FALSE),
        '.checkArgs> "fName" length error')

    expect_error(
        gG <- importNet.STRING(fName = "nonSuch.txt",
                               net = "experimental",
                               dropUnmapped = FALSE,
                               silent = TRUE,
                               writeLog = FALSE),
        '.checkArgs> "fName" error: file "nonSuch.txt" does not exist.')

    # ==== net
    expect_error(
        gG <- importNet.STRING(fName = fN,
                               net = NULL,
                               dropUnmapped = FALSE,
                               silent = TRUE,
                               writeLog = FALSE),
        '.checkArgs> "net" mode error')

    expect_error(
        gG <- importNet.STRING(fName = fN,
                               net = c("neighborhood", "fusion"),
                               dropUnmapped = FALSE,
                               silent = TRUE,
                               writeLog = FALSE),
        '.checkArgs> "net" length error')

    expect_error(
        gG <- importNet.STRING(fName = fN,
                               net = "noSuchNet",
                               dropUnmapped = FALSE,
                               silent = TRUE,
                               writeLog = FALSE),
        "Request error:")

    # ==== dropUnmapped
    expect_error(
        gG <- importNet.STRING(fName = fN,
                               net = "experimental",
                               dropUnmapped = NULL,
                               silent = TRUE,
                               writeLog = FALSE),
        '.checkArgs> "dropUnmapped" mode error')

    # ==== silent
    expect_error(
        gG <- importNet.STRING(fName = fN,
                               net = "experimental",
                               dropUnmapped = FALSE,
                               silent = NULL,
                               writeLog = FALSE),
        '.checkArgs> "silent" mode error')

    # ==== writeLog
    expect_error(
        gG <- importNet.STRING(fName = fN,
                               net = "experimental",
                               dropUnmapped = FALSE,
                               silent = TRUE,
                               writeLog = NULL),
        '.checkArgs> "writeLog" mode error')

    # ==== no logfile should have been written
    expect_false(file.exists(logName))

})


test_that("importNet.STRING() produces gG from STRING data", {
    fN <- "dataSTRING.txt"
    gG <- importNet.STRING(fName = fN,
                           net = "experimental",
                           dropUnmapped = FALSE,
                           silent = TRUE,
                           writeLog = TRUE)
    expect_equal(igraph::vcount(gG), 5)         # correct number of vertices
    expect_equal(igraph::ecount(gG), 4)         # correct number of edges

    # correct metadata
    expect_equal(attr(gG, "type"), "gG")                            # type
    expect_equal(attr(gG, "version"), "1.0")                        # version
    expect_equal(.checkArgs(attr(gG, "UUID"), "UUID"), character()) # UUID

    # correct log file
    thisLog <- readLines(logName)
    expect_true(grepl("event \\| title  \\| importNet.STRING",     thisLog[1]))
    expect_true(grepl("event \\| time",                            thisLog[2]))

    # call
    expect_true(grepl("^event \\| call   \\| importNet.STRING\\(", thisLog[3]))
    expect_true(grepl(paste0("fname = \"", fN, "\""),              thisLog[3]))
    expect_true(grepl("net = \"experimental\"",                    thisLog[3]))
    expect_true(grepl("cutoffType = \"xN\"",                       thisLog[3]))
    expect_true(grepl("val = 10000",                               thisLog[3]))
    expect_true(grepl("taxID = \"9606\"",                          thisLog[3]))
    expect_true(grepl("dropUnmapped = FALSE",                      thisLog[3]))
    expect_true(grepl("silent = TRUE",                             thisLog[3]))
    expect_true(grepl("writeLog = TRUE",                           thisLog[3]))
    expect_true(grepl(")$",                                        thisLog[3]))

    expect_true(grepl("event | note   | Read 4 edges from file.",  thisLog[4]))
    expect_true(grepl("Selected 4 edges via cutoff.",              thisLog[5]))
    expect_true(grepl("gG object has 5 vertices and 4 edges.",     thisLog[6]))
    expect_true(grepl("\"gG\"",                                    thisLog[7]))
    expect_true(all(grepl("event \\| output \\| attribute \\|", thisLog[8:11])))
    expect_true(grepl("class        | \"igraph\"",                 thisLog[8]))
    expect_true(grepl("type         | \"gG\"",                     thisLog[9]))
    expect_true(grepl("version      | \"1.0\"",                    thisLog[10]))
    expect_true(grepl("UUID",                                      thisLog[11]))
    expect_true(grepl("event | end",                               thisLog[12]))
    expect_true(grepl("^$",                                        thisLog[13]))
})


test_that("importNet.STRING() doesn't leak output if silent = TRUE", {
    fN <- "dataSTRING.txt"
    testF <- tempfile()
    capture.output(gG <- importNet.STRING(fName = fN,
                                          net = "experimental",
                                          dropUnmapped = FALSE,
                                          silent = TRUE,
                                          writeLog = FALSE),
                   file = testF)
    expect_equal(length(readLines(testF)), 0)
    expect_true(file.remove(testF))
})


test_that("importNet.STRING() doesn't produce gG from corrupt data", {
    # ToDo:
    #  - stop on missing NL on last line
    #  - stop on NA in ID
    #  - stop on NA in requested network
    #  - stop on mismatch between header elements and data columns
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
