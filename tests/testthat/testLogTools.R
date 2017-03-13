# testLogTools.R
#
#
context("Tools for reading and writing log files")

# Save the original ...
OLOG <- as.character(getOption("rete.logfile"))
testPath <- tempdir()

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

testFName <- paste("rete_", Sys.Date(), ".1.log", sep = "")
test_that("logFileName works as expected with reasonable arguments", {
    # test missing path and filename
    expect_equal(logFileName(), file.path(getwd(), testFName))
    # test missing path
    expect_equal(logFileName(fName = "x.log"), file.path(getwd(), "x.log"))
    # test missing filename
    expect_equal(logFileName(fPath = testPath), file.path(testPath, testFName))
    # test empty string filename
    expect_equal(logFileName(fName = ""), file.path(getwd(), testFName))
    # test removal of trailing /
    expect_equal(logFileName(fPath = paste(getwd(), "/", sep = "")),
                 file.path(getwd(), testFName))
    # test increment of logfilename
    writeLines("test\n", file.path(testPath, testFName))
    expect_equal(logFileName(fPath = testPath),
                 file.path(testPath, gsub("\\.1\\.log", ".2.log", testFName)))
    # cleanup, but might as well test ...
    expect_true(file.remove(file.path(testPath, testFName)))
    # test update of rete.logfile
    options("rete.logfile" = "no/such.log")
    logFileName(fPath = testPath, setOption = TRUE)
    expect_true(unlist(getOption("rete.logfile")) ==
                 file.path(testPath, testFName))
})

test_that("logMessage() rejects erroneous arguments", {
    # set rete.logfile
    logFileName(fPath = testPath, setOption = TRUE)
    # confirm that rete.logfile does not yet exist when we enter the test
    expect_false(file.exists(unlist(getOption("rete.logfile"))))
    # test missing message
    expect_error(logMessage())
    # test NULL message
    expect_error(logMessage(NULL))
    # test non-character message
    expect_error(logMessage(TRUE))
    expect_error(logMessage(as.matrix(c("a", "b"))))
    # confirm that rete.logfile was not created during the failed tests
    expect_false(file.exists(unlist(getOption("rete.logfile"))))
})

test_that("logMessage works as expected with reasonable arguments", {
    fn <- unlist(getOption("rete.logfile"))
    # confirm that rete.logfile does not yet exist when we enter the test
    expect_false(file.exists(fn))
    # test one message, file creation, correct addition of \n if \n is missing
    logMessage("a")
    expect_true(file.exists(fn))  # log file created
    #expect_equal(readChar(fn, file.info(fn)$size), "a\n")
    # test two message elements
    logMessage(c("b", "c"))
    #expect_equal(readChar(fn, file.info(fn)$size), "a\nb\nc\n")
    # test no addition of \n if \n is already there
    logMessage("d\n")
    #expect_equal(readChar(fn, file.info(fn)$size), "a\nb\nc\nd\n")
    # test logfile can be removed
    # cleanup, but might as well test ...
    expect_true(file.remove(fn))
})




# Cleanup after testing
options("rete.logfile" = OLOG)

# [END]
