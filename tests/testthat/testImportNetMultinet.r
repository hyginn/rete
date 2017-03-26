# testImportNetMultinet.R

# This is the unit test for importNet.Multinet()

context("import networks")

# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()
# ==== END SETUP AND PREPARE ===================================================


# ==== importNet.Multinet() ======================================================

test_that("importNet.Multinet() rejects erroneous parameters", {
    fN <- "Multinet.interactions.network_presence.txt"
    
    # ==== fName
    expect_error(
      gG <- importNet.Multinet(fName = NULL,
                               net = "PPI",
                               verbose = TRUE,
                               silent = TRUE,
                               writeLog = FALSE),
      '.checkArgs> "fName" mode error')
    
    expect_error(
      gG <- importNet.Multinet(fName = c("fee", "fie", "fo", "fum"),
                               net = "PPI",
                               verbose = TRUE,
                               silent = TRUE,
                               writeLog = FALSE),
      '.checkArgs> "fName" length error')
    
    expect_error(
      gG <- importNet.Multinet(fName = "nonSuch.txt",
                               net = "PPI",
                               verbose = TRUE,
                               silent = TRUE,
                               writeLog = FALSE),
      '.checkArgs> "fName" error: file "nonSuch.txt" does not exist.')
    
    # ==== net
    expect_error(
      gG <- importNet.Multinet(fName = fN,
                               net = NULL,
                               verbose = TRUE,
                               silent = TRUE,
                               writeLog = FALSE),
      '.checkArgs> "net" mode error')
    
    expect_error(
      gG <- importNet.Multinet(fName = fN,
                               net = c("PPI", "GENETIC"),
                               verbose = TRUE,
                               silent = TRUE,
                               writeLog = FALSE),
      '.checkArgs> "net" length error')
    
    expect_error(
      gG <- importNet.Multinet(fName = fN,
                               net = "noSuchNet",
                               verbose = TRUE,
                               silent = TRUE,
                               writeLog = FALSE),
      "Request error:")
    
    # ==== verbose
    expect_error(
      gG <- importNet.Multinet(fName = fN,
                               net = "PPI",
                               verbose = NULL,
                               silent = TRUE,
                               writeLog = FALSE),
      '.checkArgs> "verbose" mode error')
    
    # ==== silent
    expect_error(
      gG <- importNet.Multinet(fName = fN,
                               net = "PPI",
                               verbose = TRUE,
                               silent = NULL,
                               writeLog = FALSE),
      '.checkArgs> "silent" mode error')
    
    # ==== writeLog
    expect_error(
      gG <- importNet.Multinet(fName = fN,
                               net = "PPI",
                               verbose = TRUE,
                               silent = TRUE,
                               writeLog = NULL),
      '.checkArgs> "writeLog" mode error')
    
    # ==== no logfile should have been written
    expect_false(file.exists(logName))
    
})


test_that("importNet.Multinet() produces gG from Multinet data", {
    fN <- "Multinet.interactions.network_presence.txt"
    gG <- importNet.Multinet(fName = fN,
                             net = "PPI",
                             verbose = TRUE,
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
    expect_true(grepl("event \\| title  \\| importNet.Multinet",   thisLog[1]))
    expect_true(grepl("event \\| time",                            thisLog[2]))
    
    # call
    expect_true(grepl("^event \\| call   \\| importNet.Multinet\\(", thisLog[3]))
    expect_true(grepl(paste0("fname = \"", fN, "\""),                thisLog[3]))
    expect_true(grepl("net = \"PPI\"",                               thisLog[3]))
    expect_true(grepl("verbose = True",                              thisLog[3]))
    expect_true(grepl("silent = TRUE",                               thisLog[3]))
    expect_true(grepl("writeLog = TRUE",                             thisLog[3]))
    expect_true(grepl(")$",                                          thisLog[3]))
    
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


test_that("importNet.Multinet() doesn't leak output if silent = TRUE", {
    fN <- "Multinet.interactions.network_presence.txt"
    testF <- tempfile()
    capture.output(gG <- importNet.Multinet(fName = fN,
                                            net = "PPI",
                                            verbose = TRUE,
                                            silent = TRUE,
                                            writeLog = FALSE),
                   file = testF)
    expect_equal(length(readLines(testF)), 0)
    expect_true(file.remove(testF))
})


test_that("importNet.Multinet() doesn't produce gG from corrupt data", {
    fN <- "Multinet.interactions.network_presence.txt"
    
    # ==== stop on missing NL on last line
    inputdata <- file(fN, open = "r")
    expect_error(
      tail(inputdata, 1) <- !c("\n"),
      '.checkArgs> stop on missing newline on last line')
    
    # ==== stop on NA in requested network
    network <- list("PPI", "METABOLIC", "GENETIC", "PHOSPHORYLATION",
                    "REGULATORY", "SIGNALING", "NUM_NETWORKS")
    expect_error(
      NotAvailable <- !("net" %in% names(network)),  # "net" input in importNet.Multinet()
      '.checkArgs> stop on not available in requested network')
    
    # ==== stop on mismatch between header elements and data columns
    inputdata <- file(fN, open = "r")
    oneLine <- readLines(inputdata, n = 1, warn = TRUE)
    skipheader <- oneLine[-1]  # skip the header line
    current <- 1
    characters <- list("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", 
                      "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
    numbers <- list("1", "2", "3", "4", "5", "6", "7", "8", "9", "0")
    while (length(oneLine) > 0) {
      expect_error(  # expect error if INTERACTION_NAME mismatch data columns
        oneLine[current][1] <- !c("_", characters, numbers),
        '.checkArgs> stop on mismatch between header elements and data columns')
    
      expect_error(  # expect error if network mismatch data columns
        !is.numeric(list(oneLine[current][2:])),
        '.checkArgs> stop on mismatch between header elements and data columns')
      current <- current + 1
    }
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
