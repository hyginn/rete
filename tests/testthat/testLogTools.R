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
    # platform appropriate conversions.
    if (.Platform$OS.type == "windows") {
        logMessage("a")
        expect_true(file.exists(fn))  # log file created
        expect_equal(readChar(fn, file.info(fn)$size), "a\r\n")
        # test two message elements
        logMessage(c("b", "c"))
        expect_equal(readChar(fn, file.info(fn)$size), "a\r\nb\r\nc\r\n")
        # test no extra \r\n if \r\n is already there
        logMessage("d\r\n")
        expect_equal(readChar(fn, file.info(fn)$size), "a\r\nb\r\nc\r\nd\r\n")
        # test replacement of wrong linebreak
        logMessage("\ne\n")
        expect_equal(readChar(fn, file.info(fn)$size),
                     "a\r\nb\r\nc\r\nd\r\n\r\ne\r\n")
    } else {
        logMessage("a")
        expect_true(file.exists(fn))  # log file created
        expect_equal(readChar(fn, file.info(fn)$size), "a\n")
        # test two message elements
        logMessage(c("b", "c"))
        expect_equal(readChar(fn, file.info(fn)$size), "a\nb\nc\n")
        # test no extra \n if \n is already there
        logMessage("d\n")
        expect_equal(readChar(fn, file.info(fn)$size), "a\nb\nc\nd\n")
        # test replacement of wrong linebreak
        logMessage("\r\ne\r\n")
        expect_equal(readChar(fn, file.info(fn)$size), "a\nb\nc\nd\n\ne\n")
    }

    # cleanup, but might as well test ...
    expect_true(file.remove(fn))
})

# tests for: attachUUID(object, overwrite=TRUE)
test_that("attachUUID rejects objects that don't exist", {
    # pass object which does not exist into function
    expect_error(attachUUID(x))

    # should reject NULL objects
    expect_error(attachUUID(NULL))
})

test_that("attachUUID does not reject objects that exist", {
    object <- "c"
    expect_error(attachUUID(object), NA)
})

test_that("attachUUID attaches UUID if the overwrite flag is FALSE AND the object does not have a UUID", {
    # generate an object which does not have a UUID set
    object <- "c"
    # call function with overwrite = FALSE
    object <- attachUUID(object, overwrite = FALSE)
    # expect the function to have attached a UUID to object
    # wt> can't actually do this because global object isn't modified by local function
    expect_false(is.null(attr(object, "UUID")))
})

test_that("attachUUID does not attach UUID if the overwrite flag is FALSE AND object already has UUID attribute", {
    # generate an object which has a UUID set
    object <- "c"
    generatedUUID <- uuid::UUIDgenerate()
    attr(object, "UUID") <- generatedUUID
    # call function with overwrite set to FALSE
    object <- attachUUID(object, overwrite = FALSE)
    # expect object's UUID to not have changed
    expect_equal(generatedUUID, attr(object, "UUID"))
})

test_that("attachUUID attaches the UUID if the overwrite flag is TRUE AND object does not have a UUID", {
    # generate an object which does not have a UUID set
    object <- "c"
    # call function with overwrite set to TRUE
    object <- attachUUID(object, overwrite = TRUE)
    # expect the function to have attached a UUID to object
    expect_false(is.null(attr(object, "UUID")))
})

test_that("if the overwrite flag is TRUE AND object has a UUID, should attach UUID", {
    # generate an object which has a UUID
    object <- "c"
    # store this UUID
    generatedUUID <- uuid::UUIDgenerate()
    attr(object, "UUID") <- generatedUUID
    # call function with overwrite set to TRUE
    object <- attachUUID(object, overwrite = TRUE)
    # expect the new object's UUID to be different from the old UUID
    expect_false(generatedUUID == attr(object, "UUID"))
})

# tests for: extractAttributes(object)
test_that("extractAttributes rejects objects that don't exist", {
    # pass reference of object which does not exist into function
    expect_error(extractAttributes(x))

    # expect NULL object to raise error
    expect_error(extractAttributes(NULL))
})

test_that("extractAttributes does not reject objects that exist", {
    # pass reference of object which exists into function
    x <- "c"
    expect_error(extractAttributes(x), NA)
})

test_that("extractAttributes return NULL for an object with no attributes", {
    # create object with no attributes
    x <- "c"
    # expect no attributes to be returned
    expect_equal(extractAttributes(x), NULL)
})

test_that("extractAttributes logs one attribute for an object with one attribute", {
    # create object with one attribute
    x <- "c"
    attr(x, "testAttribute") <- "testValue"
    # expect one attribute log text to be returned
    result <- extractAttributes(x)
    expect_equal(result, "attribute\t|\ttestAttribute\t|\ttestValue\n")
})

test_that("extractAttributes logs all attributes for an object with multiple attributes", {
    # create object with multiple attributes (more than one)
    # expect all of the attributes log text to be returned
    x <- "c"
    attr(x, "testAttribute1") <- "testValue1"
    attr(x, "testAttribute2") <- "testValue2"
    attr(x, "testAttribute3") <- "testValue3"
    # expect one attribute log text to be returned
    result1 <- "attribute\t|\ttestAttribute1\t|\ttestValue1\n"
    result2 <- "attribute\t|\ttestAttribute2\t|\ttestValue2\n"
    result3 <- "attribute\t|\ttestAttribute3\t|\ttestValue3\n"
    expect_equal(extractAttributes(x), paste(result1, result2, result3, sep=""))
})


# tests for: logEvent(eventTitle, eventCall, input = c(), output = c())

### bs> Finally, could the "event" just be a comment, or other message?
### bs> This needs a bit of thought.
### wt>
### wt> The "event" should come from a predetermined set of all possible events
### wt> that could happen in rete - it would be nice to have some sort of categorized dictionary
### wt> of events to outcomes.
### wt>
### bs>
### bs> Test that messages are separated from each other by blank lines.
### wt> Ok.

test_that("logEvent raises an error if the event title is not passed in", {
    # call function without eventTitle parameter
    # expect function to error
})

test_that("logEvent raises an error if the event call is not passed in", {
    # call function without call
    # expect function to error
})

test_that("logEvent is able to accept both objects or filenames in list of input and output parameters", {
    # call function with name that does not exist as a object or filename
    # function should error

    # call function with name that exists as an object but not as a filename
    # function should not error

    # call function with name that does not exist as an object but exists as a file
    # function should not error

    # call function with name that exists as both object or as filename
    # function should use object and should not error
})

test_that("logEvent tests that the event title comes from a predefined set of events", {

})

test_that("logEvent separates messages from each other using blank lines", {

})

test_that("logEvent should log the date and time of the event", {
    # call function with arbitrary event and no input or outputs
    # expect function to have logged the event with date and time
    # expect function to have logged the event title
    # expect function to have logged no attributes of any input or output objects
})

# only one test case for input attributes since already tested above
test_that("logEvent should log the input attributes, if any", {
    # create potential input object with two attributes
    # call function with arbitrary event
    # expect function to have logged the event with date and time
    # expect function to have logged the event title
    # expect function to have logged both input attributes
})

test_that("logEvent should log the output attributes, if any", {
    # create potential output object with two attributes
    # call function with arbitrary event
    # expect function to have logged the event with date and time
    # expect function to have logged the event title
    # expect function to have logged both output attributes
})

test_that("logEvent should log both input and output attributes", {
    # create potential input object with two attributes
    # create potential output object with two attributes
    # call function with arbitrary event
    # expect function to have logged the event with date and time
    # expect function to have logged the event title
    # expect function to have logged both input attributes
    # expect function to have logged both output attributes
})


# tests for: findUUID(uuid, path)
test_that("function errors if an invalid path has been provided", {
    # call function with invalid file path
    # expect function to error
})

test_that("if the uuid does not occur in the logs, that nothing is returned", {
    # set up log file with uuids that will not be called by the function test
    # call function with uuid which does not occur in the logs
    # expect function to not return any instances
    ### bs> What should the function do instead? Return NULL, NA, or ""?
    ### wt> The function should instead return NULL.
})

test_that("one instance is returned if uuid occurs in the logs once", {
    # set up log file with one instance of a uuid which will be called by the function
    # call function with same uuid
    # expect function to return one instance
    ### bs> What _is_ actually returned? Just the line? The entire event block?
    ### bs>    The filename too?
})

test_that("multiple instances are return if uuid occurs in the logs more than once", {
    # set up log file with 3 instances of a uuid which will be called by the function
    # call function with same uuid
    # expect function to return the 3 instances of uuid being mentioned, in chronological order
})

###

# tests for: getProvenance(object, path)

### bs> I think you should pass an UUID as an argument, not an object. Yes?
### wt> Yes. Will change pseudocode.

test_that("function throws error if object does not exist", {
    # call function with object which does not exist
    # expect function to error
    ### bs> IF you pass an UUID, your function needs to test whether the object is structured like an UUID.
    ### bs> you should use the .checkArgs() function for that and I will add a
    ### bs> keyword directive    like = "UUID"  to .checkArgs() . Good?
    ### wt> great, thank you.
})

test_that("function throws error if object's uuid does not exist in the logs", {
    # set up mock log file with uuids that will not be called by the function test
    # call function with object uuid which does not occur in the logs
    # expect function to error
    ### bs> I would not call this an error. I think this test is actually
    ### bs> the same as the second test in this group.
    ### wt> I will remove this test.
})

test_that("function throws error if a path specified is not a directory", {
    # call function with invalid path
    # expect function to error
    ### bs> yes. .checkArgs has an option for that
})

### bs> I think it would be useful to also throw an error if there are no log files at all in
### bs> the directory since that's almost certainly not intended.
### wt> yes, will add pseudocode.

test_that("given an object whose UUID occurs in the logs, returns the correct provenance", {
    # set up multiple log events where a particular object was manipulated, from a filename without a UUID
    # call function with valid object and directory
    # expect function to return the intermediaries of UUIDs and events which lead up to the current state of the object
    ### bs> Here too: you need to decide what is being returned. Probably full event blocks plus filenames.
})

### bs> provenance can be branching since objects can have multiple parents. Need to test that these are correctly queued
### bs> and handled. Suggest depth-first, not breadth-first.

# Cleanup after testing
options("rete.logfile" = OLOG)

# [END]
