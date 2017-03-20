# testLogTools.R
#
#
context("Tools for reading and writing log files")

# Save the original ...
OLOG <- as.character(getOption("rete.logfile"))
NL <- .PlatformLineBreak()


# ==== logFileName() ===========================================================

testPath <- tempdir()
testFName <- paste("rete_", Sys.Date(), ".1.log", sep = "")

test_that("logFileName() rejects erroneous arguments", {
    # test non-existing path
    expect_error(logFileName(fPath = "no/such/path"))
    # test NULL path
    expect_error(logFileName(fPath = NULL))
    # test null-length path
    expect_error(logFileName(fPath = character()))
    # test vector of paths
    expect_error(logFileName(fPath = c("~", "~")))

    # test NULL filename
    expect_error(logFileName(fPath = testPath, fName = NULL))
    # test null-length filename
    expect_error(logFileName(fPath = testPath, fName = character()))
    # test vector of filenames
    expect_error(logFileName(fPath = testPath, fName = c("a", "b")))

    # test option of wrong type
    expect_error(logFileName(setOption = 1))
    # test NULL option
    expect_error(logFileName(setOption = NULL))
    # test null-length option
    expect_error(logFileName(setOption = logical()))
    # test vector of options
    expect_error(logFileName(setOption = c(TRUE, TRUE)))
})


test_that("logFileName works as expected with reasonable arguments", {

    # test missing path becomes getwd() or getwd()/logs
    if (dir.exists(file.path(getwd(), "logs"))) {
        expect_equal(dirname(logFileName()), file.path(getwd(), "logs"))
    } else {
        expect_equal(dirname(logFileName()), getwd())
    }

    # remove ".log" files from testPath if any exist
    if (length(list.files(path = testPath,
                          all.files = TRUE,
                          pattern = "\\.log$")) > 0 ) {
        invisible(file.remove(list.files(path = testPath,
                                         all.files = TRUE,
                                         pattern = "\\.log$",
                                         full.names = TRUE)))
    }

    # initial state: no ".log" files in testPath
    expect_true(length(list.files(path = testPath,
                                  all.files = TRUE,
                                  pattern = "\\.log$")) == 0)

    # Actual tests:
    # test missing path becomes getwd()
    expect_equal(logFileName(fName = "x.log"), file.path(getwd(), "x.log"))
    # test missing filename
    expect_equal(logFileName(fPath = testPath), file.path(testPath, testFName))
    # test empty string filename
    expect_equal(logFileName(fName = ""), file.path(getwd(), testFName))
    # test removal of trailing / from path
    expect_equal(logFileName(fPath = paste(getwd(), "/", sep = "")),
                 file.path(getwd(), testFName))
    # test increment of logfilename
    # create a file with a log-file pattern
    writeLines("test", file.path(testPath, testFName))
    expect_equal(logFileName(fPath = testPath),
                 file.path(testPath, gsub("\\.1\\.log", ".2.log", testFName)))
    # cleanup, but might as well test ...
    expect_true(file.remove(file.path(testPath, testFName)))

    # test update of rete.logfile
    options("rete.logfile" = "no/such.log") # break global filename ..
    invisible(logFileName(fPath = testPath, setOption = TRUE)) # fix it ..
    expect_true(unlist(getOption("rete.logfile")) ==
                 file.path(testPath, testFName))

    # final state:: no ".log" files in testPath
    expect_true(length(list.files(path = testPath,
                                  all.files = TRUE,
                                  pattern = "\\.log$")) == 0)

})


# ==== logMessage() tests ======================================================

# set up fresh state
testPath <- tempdir()
testFName <- paste("rete_", Sys.Date(), ".1.log", sep = "")
if (dir.exists(file.path(testPath, "logs"))) {
    if (file.exists(logFileName(fPath = testPath))) {
        file.remove(logFileName(fPath = testPath))
    }
    file.remove(file.path(testPath, "logs"))
}

if (file.exists(logFileName(fPath = testPath))) {
    file.remove(logFileName(fPath = testPath))
}

invisible(logFileName(fPath = testPath, setOption = TRUE))


test_that("logMessage() rejects erroneous arguments", {
    # get the log file name
    FN <- unlist(getOption("rete.logfile"))

    # confirm that rete.logfile does not yet exist when we enter the test
    expect_false(file.exists(FN))
    # test missing message
    expect_error(logMessage())
    # test NULL message
    expect_error(logMessage(NULL))
    # test non-character message
    expect_error(logMessage(TRUE))
    expect_error(logMessage(as.matrix(c("a", "b"))))
    # confirm that rete.logfile was not created during the failed tests
    expect_false(file.exists(FN))
})


test_that("logMessage works as expected with reasonable arguments", {
    # get the log file name
    FN <- unlist(getOption("rete.logfile"))

    # confirm that rete.logfile does not yet exist when we enter the test
    expect_false(file.exists(FN))

    # test one message, file creation, correct addition of \n if \n is missing
    logMessage("a")
    expect_true(file.exists(FN))  # log file created
    expect_equal(readChar(FN, file.info(FN)$size),
                 paste0("a", NL))
    # test two message elements
    logMessage(c("b", "c"))
    expect_equal(readChar(FN, file.info(FN)$size),
                 paste0("a", NL, "b", NL, "c", NL))
    # test no extra \n if \n is already there
    logMessage(paste0("d", NL))
    expect_equal(readChar(FN, file.info(FN)$size),
                 paste0("a", NL, "b", NL, "c", NL, "d", NL))
    # test correct empty line
    logMessage("")
    expect_equal(readChar(FN, file.info(FN)$size),
                 paste0("a", NL, "b", NL, "c", NL, "d", NL, NL))
    # test replacement of wrong internal linebreak
    logMessage("e\nf\r\ng\n")
    expect_equal(readChar(FN, file.info(FN)$size),
                 paste0("a", NL, "b", NL, "c", NL, "d", NL, NL,
                        "e", NL, "f", NL, "g", NL))

    # cleanup, but might as well test ...
    expect_true(file.remove(FN))
})


# ==== getUUID() ===============================================================

# pattern to grep for UUIDs
patt<-"[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}"

test_that("getUUID rejects inappropriate objects", {

    expect_error(getUUID("testAsset_01"),
                 'Object "testAsset_01" does not exist.')
    "testAsset_01" <- NULL
    expect_error(getUUID("testAsset_01"),
                 'Object "testAsset_01" is NULL.')
    rm("testAsset_01")
})

test_that("getUUID gets an UUID for an object that has none", {
    testAsset_02 <- "a"
    myUUID <- getUUID("testAsset_02")
    expect_true(grepl(patt, myUUID))
    rm("testAsset_02")
})

test_that("getUUID returns the old UUID if the object has one", {
    testAsset_03 <- "b"
    attr(testAsset_03, "UUID") <- "fakeUUID"
    expect_equal(getUUID("testAsset_03"), "fakeUUID")
    rm("testAsset_03")
})

test_that("getUUID returns a new UUID if overwrite is TRUE", {
    testAsset_04 <- "c"
    attr(testAsset_04, "UUID") <- "fakeUUID"
    expect_true(grepl(patt, getUUID("testAsset_04", overwrite = TRUE)))
    rm("testAsset_04")
})

# ==== .extractAttributes() ====================================================

if (exists("tmp")) { rm(tmp) }
test_that(".extractAttributes() rejects an object that does not exist", {
    expect_error(.extractAttributes("tmp", role = "input"),
                 "Object \"tmp\" does not exist.")
})

test_that(".extractAttributes() rejects NULL objects", {
    tmp <- NULL
    expect_error(.extractAttributes("tmp", role = "input"),
                 "Object \"tmp\" is NULL.")
})

test_that(".extractAttributes() returns character() if object has none", {
    tmp <- "c"
    expect_equal(.extractAttributes("tmp", role ="input"), character())
})

test_that(".extractAttributes() throws error for missing or incorrect role", {
    tmp <- "c"
    expect_error(.extractAttributes("tmp"),
                 "role parameter must be provided.")
    expect_error(.extractAttributes("tmp", role = "nonsuch"))
})

test_that(".extractAttributes() logs one attribute", {
    tmp <- "c"
    attr(tmp, "x") <- "u"
    expect_equal(.extractAttributes("tmp", role = "input"),
                 "event | input  | attribute | x            | \"u\"")
    attr(tmp, "x") <- c("u", "v")
    expect_equal(.extractAttributes("tmp", role = "output"),
                 "event | output | attribute | x            | (u, v)")
    attr(tmp, "x") <- LETTERS
    expect_equal(.extractAttributes("tmp", role = "input"),
        "event | input  | attribute | x            | (A, B, C, ... (26) )")
})

test_that(".extractAttributes() logs multiple attributes", {
    tmp <- "c"
    attr(tmp, "x") <- "u"
    attr(tmp, "y") <- c("u", "v")
    result <- .extractAttributes("tmp", role = "input")

    expect_equal(result[1],
                 "event | input  | attribute | x            | \"u\"")
    expect_equal(result[2],
                 "event | input  | attribute | y            | (u, v)")

    result <- .extractAttributes("tmp", role = "output")
    expect_equal(result[1],
                 "event | output | attribute | x            | \"u\"")
    expect_equal(result[2],
                 "event | output | attribute | y            | (u, v)")
})


# ==== logEvent() ==============================================================

logFileName(testPath, testFName, setOption = TRUE)
FN <- unlist(getOption("rete.logfile"))
if (file.exists(FN)) { file.remove(FN)}

test_that("logEvent raises an error if title and/or call are missing", {
    expect_error(logEvent())
    expect_error(logEvent(eventTitle = "test"))
    expect_error(logEvent(eventCall = "f1(p = P)"))
})

test_that("logEvent raises an error if any parameter is not mode character", {
    expect_error(logEvent(eventTitle = TRUE,
                          eventCall = "f1(p = P)",
                          input = "this",
                          notes = "note",
                          output = "that"))
    expect_error(logEvent(eventTitle = "test",
                          eventCall = TRUE,
                          input = "this",
                          notes = "note",
                          output = "that"))
    expect_error(logEvent(eventTitle = "test",
                          eventCall = "f1(p = P)",
                          input = TRUE,
                          notes = "note",
                          output = "that"))
    expect_error(logEvent(eventTitle = "test",
                          eventCall = "f1(p = P)",
                          input = "this",
                          notes = TRUE,
                          output = "that"))
    expect_error(logEvent(eventTitle = "test",
                          eventCall = "f1(p = P)",
                          input = "this",
                          notes = "note",
                          output = TRUE))
})

test_that("logEvent raises an error if objects don't exist or are NULL", {
    expect_error(logEvent(eventTitle = "test",
                          eventCall = "f1(p = P)",
                          input = "this"))
    expect_error(logEvent(eventTitle = "test",
                          eventCall = "f1(p = P)",
                          output = "that"))
    this <- "something"
    expect_error(logEvent(eventTitle = "test",
                          eventCall = "f1(p = P)",
                          input = c("this", "that")))
    expect_error(logEvent(eventTitle = "test",
                          eventCall = "f1(p = P)",
                          output = c("this", "that")))
    that <- NULL
    expect_error(logEvent(eventTitle = "test",
                          eventCall = "f1(p = P)",
                          input = c("this", "that")))
    expect_error(logEvent(eventTitle = "test",
                          eventCall = "f1(p = P)",
                          output = c("this", "that")))
})


test_that("The previous tests did not lead to creating a logfile", {
    expect_false(file.exists(FN))
})


test_that("logEvent writes a structured message", {

    FN <- unlist(getOption("rete.logfile"))
    if (file.exists(FN)) { file.remove(FN)}

    testAsset_04 <- TRUE
    attr(testAsset_04, "UUID") <- uuid::UUIDgenerate()

    testAsset_05 <- 1:10
    attr(testAsset_05, "att") <- "a1"
    attr(testAsset_05, "UUID") <- uuid::UUIDgenerate()

    testAsset_06 <- letters
    names(testAsset_06) <- LETTERS
    attr(testAsset_06, "UUID") <- uuid::UUIDgenerate()

    logEvent(eventTitle = "test 1",
             eventCall = "f1(p = P)",
             input = c("testAsset_04", "testAsset_05"),
             notes = c("lorem ipsum", "dolor sit amet"),
             output = "testAsset_06")

    expect_true(file.exists(FN))

    logResult <- readLines(FN)

    expect_equal(logResult[1],  "event | title  | test 1")
    expect_true(grepl("event | time   | [0-9 -:]$",
                      logResult[2]))
    expect_equal(logResult[3],  "event | call   | f1(p = P)")
    expect_equal(logResult[4],  "event | input  | \"testAsset_04\"")
    expect_true(grepl("event | input  | attribute | UUID", logResult[5]))
    expect_true(grepl(patt, logResult[5]))
    expect_equal(logResult[6],  "event | input  | \"testAsset_05\"")
    expect_equal(logResult[7],
                 "event | input  | attribute | att          | \"a1\"")
    expect_true(grepl("event | input  | attribute | UUID", logResult[8]))
    expect_true(grepl(patt, logResult[8]))
    expect_equal(logResult[9],  "event | note   | lorem ipsum")
    expect_equal(logResult[10], "event | note   | dolor sit amet")
    expect_equal(logResult[11], "event | output | \"testAsset_06\"")
    expect_equal(logResult[12],
        "event | output | attribute | names        | (A, B, C, ... (26) )")
    expect_true(grepl("event | output | attribute | UUID", logResult[13]))
    expect_true(grepl(patt, logResult[13]))
    expect_equal(logResult[14], "event | end")
    expect_equal(logResult[15], "")

    # append the same again
    logEvent(eventTitle = "test 2",
             eventCall = "f1(p = P)",
             input = c("testAsset_04", "testAsset_05"),
             notes = c("lorem ipsum", "dolor sit amet"),
             output = "testAsset_06")

    logResult <- readLines(FN)
    expect_equal(length(logResult), 30)

    rm("testAsset_04")
    rm("testAsset_05")
    rm("testAsset_06")
    file.remove(FN)
})

# ==== findUUID() ==============================================================


test_that("findUUID() rejects invalid arguments", {
    tmpUUID <- uuid::UUIDgenerate()
    expect_error(findUUID())
    expect_error(findUUID(uuid = tmpUUID, logDir = "no/such/path"))
    expect_error(findUUID(uuid = tmpUUID, logDir = NULL))
    expect_error(findUUID(uuid = tmpUUID, logDir = character()))
    expect_error(findUUID(uuid = tmpUUID, logDir = c("~", "~")))
    expect_error(findUUID(uuid = "not a UUID"))
    expect_error(findUUID(uuid = tmpUUID, ext = NULL))
    expect_error(findUUID(uuid = tmpUUID, recursive = NULL))
})

test_that("findUUID() stops if no log files in logDir", {
    expect_error(findUUID(uuid = tmpUUID, ext = "no.such.log$"))
})

test_that("findUUID() returns the correct number of events", {

    FN <- unlist(getOption("rete.logfile"))
    if (file.exists(FN)) { file.remove(FN)}

    testAsset_04 <- TRUE
    attr(testAsset_04, "UUID") <- uuid::UUIDgenerate()

    testAsset_05 <- 1:10
    attr(testAsset_05, "att") <- "a1"
    attr(testAsset_05, "UUID") <- uuid::UUIDgenerate()

    testAsset_06 <- letters
    names(testAsset_06) <- LETTERS
    attr(testAsset_06, "UUID") <- uuid::UUIDgenerate()

    logEvent(eventTitle = "test 1",
             eventCall = "doAThing(thing)",
             input = c("testAsset_04", "testAsset_05"),
             notes = c("lorem ipsum", "dolor sit amet"),
             output = "testAsset_06")

    logEvent(eventTitle = "test 2",
             eventCall = "doAnotherThing(nextThing)",
             input = c("testAsset_04", "testAsset_05"),
             notes = c("lorem ipsum", "dolor sit amet"),
             output = "testAsset_06")

    a4UUIDold <- attr(testAsset_04, "UUID")
    a4UUIDnew <- uuid::UUIDgenerate()
    attr(testAsset_04, "UUID") <- a4UUIDnew

    a5UUID    <- attr(testAsset_05, "UUID")
    xxUUID    <- uuid::UUIDgenerate()

    logEvent(eventTitle = "test 3",
             eventCall = "doALastThing(finalThing)",
             input = c("testAsset_04", "testAsset_05"),
             output = "testAsset_06")

    # the event log
    #  - contains no mentions of xxUUID
    #  - one mention of a4UUIDnew
    #  - two mentions of a4UUIDold
    #  - three mentions of a5UUID

    expect_equal(length(findUUID(uuid = xxUUID)), 0)
    expect_equal(length(findUUID(uuid = a4UUIDnew)), 1)
    expect_equal(length(findUUID(uuid = a4UUIDold)), 2)
    expect_equal(length(findUUID(uuid = a5UUID)), 3)

    # ...that the events have the correct boundaries
    events <- findUUID(uuid = a5UUID)

    expect_true(grepl("^event | title  | test", events[1]))
    expect_true(grepl("^event | end$", events[1]))
    expect_true(grepl("^event | title  | test", events[2]))
    expect_true(grepl("^event | end$", events[2]))
    expect_true(grepl("^event | title  | test", events[3]))
    expect_true(grepl("^event | end$", events[3]))

    rm("testAsset_04")
    rm("testAsset_05")
    rm("testAsset_06")
    file.remove(FN)
})


# ==== getProvenance() =========================================================


#
# # tests for: getProvenance(uuid, uuidPath)
#
# test_that("function throws error if object does not exist", {
#     # call function with object which does not exist
#     # expect function to error
#     ### bs> IF you pass an UUID, your function needs to test whether the
#              object is structured like an UUID.
#     ### bs> you should use the .checkArgs() function for that and I will add a
#     ### bs> keyword directive    like = "UUID"  to .checkArgs() . Good?
#     ### wt> great, thank you.
# })
#
# test_that("function throws error if object's uuid does not exist
#            in the logs", {
#     # set up mock log file with uuids that will not be called by the
#       function test
#     # call function with object uuid which does not occur in the logs
#     # expect function to error
#     ### bs> I would not call this an error. I think this test is actually
#     ### bs> the same as the second test in this group.
#     ### wt> I will remove this test.
# })
#
# test_that("function throws error if a path specified is not a directory", {
#     # call function with invalid path
#     # expect function to error
#     ### bs> yes. .checkArgs has an option for that
# })
#
# ### bs> I think it would be useful to also throw an error if there are no log
#         files at all in
# ### bs> the directory since that's almost certainly not intended.
# ### wt> yes, will add pseudocode.
#
# test_that("given an object whose UUID occurs in the logs, returns the
#           correct provenance", {
#     # set up multiple log events where a particular object was manipulated,
#      from a filename without a UUID
#     # call function with valid object and directory
#     # expect function to return the intermediaries of UUIDs and events which
#      lead up to the current state of the object
#     ### bs> Here too: you need to decide what is being returned. Probably full
#      event blocks plus filenames.
# })
#
# ### bs> provenance can be branching since objects can have multiple parents.
#  Need to test that these are correctly queued
# ### bs> and handled. Suggest depth-first, not breadth-first.

# Cleanup after testing
FN <- unlist(getOption("rete.logfile"))
if (file.exists(FN)) { file.remove(FN)}

options("rete.logfile" = OLOG)

# [END]
