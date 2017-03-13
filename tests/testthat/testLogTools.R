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
    expect_error(extractAttributes(x), "object 'x' not found")

    # expect NULL object to raise error
    expect_error(extractAttributes(NULL), "object is null")
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
    expect_equal(result, "event\t|\tinput\t|\tattribute\t|\ttestAttribute\t|\ttestValue\n")
})

test_that("extractAttributes logs all attributes for an object with multiple attributes", {
    # create object with multiple attributes (more than one)
    # expect all of the attributes log text to be returned
    x <- "c"
    attr(x, "testAttribute1") <- "testValue1"
    attr(x, "testAttribute2") <- "testValue2"
    attr(x, "testAttribute3") <- "testValue3"
    # expect one attribute log text to be returned
    result1 <- "event\t|\tinput\t|\tattribute\t|\ttestAttribute1\t|\ttestValue1\n"
    result2 <- "event\t|\tinput\t|\tattribute\t|\ttestAttribute2\t|\ttestValue2\n"
    result3 <- "event\t|\tinput\t|\tattribute\t|\ttestAttribute3\t|\ttestValue3\n"
    expect_equal(extractAttributes(x), paste(result1, result2, result3, sep=""))
})


# tests for: logEvent(eventTitle, eventCall, input = list(), output = list(), fPath = getwd())

### bs> Finally, could the "event" just be a comment, or other message?
### bs> This needs a bit of thought.
### wt>
### wt> The "event" should come from a predetermined set of all possible events
### wt> that could happen in rete - it would be nice to have some sort of categorized dictionary
### wt> of events to outcomes.

test_that("logEvent raises an error if the event title is not passed in", {
    # call function without eventTitle parameter
    eventCall <- "function1(param1, param2)"
    # expect function to error
    expect_error(logEvent(eventCall = eventCall), "eventTitle is NULL")
})

test_that("logEvent raises an error if the event call is not passed in", {
    # call function without call
    eventTitle <-"Graph generation"
    # expect function to error
    expect_error(logEvent(eventTitle = eventTitle), "eventCall is NULL")
})

test_that("logEvent raises an error if both the event call and event title are not passed in", {
    expect_error(logEvent())
})

test_that("logEvent is able to accept both objects or filenames in list of input and output parameters", {
    logFile <- unlist(getOption("rete.logfile"))
    expect_false(file.exists(logFile))

    eventTitle <-"Graph generation"
    eventCall <- "function1(param1, param2)"
    # call function with input name that DNE as an object or filename
    # function should error
    expect_error(logEvent(eventTitle, eventCall, input = list(x)), "object 'x' not found")

    # call function with output name that DNE as object or filename
    expect_error(logEvent(eventTitle, eventCall, output = list(x)), "object 'x' not found")

    # call function with NULL input
    # function should error
    expect_error(logEvent(eventTitle, eventCall, input = list(NULL)), "input contains NULL object")

    # call function with NULL output
    # function should error
    expect_error(logEvent(eventTitle, eventCall, output = list(NULL)), "output contains NULL object")

    # call function with input name that exists as an object but not as a filename
    # function should not error
    testInputObject1 <- c("testparam1", "testparam2", "testparam3")
    attr(testInputObject1, "testInputAttr1") <- "testInputValue1"
    expect_error(logEvent(eventTitle, eventCall, input = list(testInputObject1)), NA)

    # call function with output name that exists as an object but not as a filename
    testOutputObject1 <- c("testparam1", "testparam2", "testparam3")
    attr(testOutputObject1, "testOutputAttr1") <- "testOutputValue1"
    expect_error(logEvent(eventTitle, eventCall, output = list(testOutputObject1)), NA)

    # call function with input name that does not exist as an object but exists as a file
    # function should not error
    testInputFile1 <- "testInputFile1.rds"
    save(testInputFile1, file = "testInputFile1.rds")
    expect_error(logEvent(eventTitle, eventCall, input = list(testInputFile1)), NA)
    expect_true(file.remove(testInputFile1))

    # call function with output name that DNE as an object but exists as a file
    # function should not error
    testOutputFile1 <- "testOutputFile1.rds"
    save(testOutputFile1, file = "testOutputFile1.rds")
    expect_error(logEvent(eventTitle, eventCall, output = list(testOutputFile1)), NA)
    expect_true(file.remove(testOutputFile1))

    expect_true(file.remove(logFile))

})

test_that("logEvent tests that the event title comes from a predefined set of events", {
    eventTitle <-"Graph generation"
    eventCall <- "function1(param1, param2)"
})

test_that("logEvent writes proper messages", {
    logFile <- unlist(getOption("rete.logfile"))
    expect_false(file.exists(logFile))

    eventTitle <-"Graph generation"
    eventCall <- "function1(param1, param2)"

    testInputObject1 <- c("testparam1", "testparam2", "testparam3")
    attr(testInputObject1, "testBlankAttr1") <- "testBlankValue1"

    testInputObject2 <- c("testparam4")
    attr(testInputObject2, "testBlankAttr2") <- "testBlankValue2"
    expect_error(logEvent(eventTitle, eventCall, input = list(testInputObject1, testInputObject2)), NA)

    expect_true(file.exists(logFile))

    result1 <- "event\t|\tinput\t|\tattribute\t|\ttestBlankAttr1\t|\ttestBlankValue1\n"
    result2 <- "event\t|\tinput\t|\tattribute\t|\ttestBlankAttr2\t|\ttestBlankValue2\n"
    call <- "event\t|\tcall\t|\tfunction1(param1, param2)\n"
    title <- "comment\t|\ttitle\t|\tGraph generation\n"
    dateTime <- paste("comment\t|\tdateTime\t|\t", Sys.time(), "\n", sep = "")

    expect_equal(readChar(logFile, file.info(logFile)$size), paste(result1, result2, call, title, dateTime, sep = ""))
    expect_true(file.remove(logFile))
})


# tests for: findUUID(uuid, uuidPath)
test_that("findUUID errors if an invalid path has been provided", {
    randomUUID <- "34cb8d64-4422-4232-9d58-8521b5a90ada"

    # Professor Steipe's tests from testLogTools.R
    # test non-existent path
    expect_error(findUUID(randomUUID, uuidPath = "no/such/path"))
    # test NULL path
    expect_error(findUUID(randomUUID, uuidPath = NULL))
    # test null-length path
    expect_error(findUUID(randomUUID, uuidPath = character()))
    # test vector of paths
    expect_error(findUUID(randomUUID, uuidPath = c("~", "~")))
})

test_that("findUUID errors if UUID is not valid", {
    # test random string
    badUUID1 <- "abc"
    expect_error(findUUID(badUUID1))
    # test vector
    badUUID2 <- c("not", "good")
    expect_error(findUUID(badUUID2))
    # test something very close
    badUUID3 <- "89337ch-389-3938-9383-938e3cu0938"
    expect_error(findUUID(badUUID3))
})

test_that("if the uuid does not occur in the logs, that nothing is returned", {
    # set up log file with uuids that will not be called by the function test
    logFile <- unlist(getOption("rete.logfile"))
    expect_false(file.exists(logFile))

    eventTitle <-"Graph generation"
    eventCall <- "function1(param1, param2)"

    testInputObject1 <- c("testparam1", "testparam2", "testparam3")
    attr(testInputObject1, "testBlankAttr1") <- "testBlankValue1"
    testInputObject1 <- attachUUID(testInputObject1)

    testInputObject2 <- c("testparam4")
    attr(testInputObject2, "testBlankAttr2") <- "testBlankValue2"
    testInputObject2 <- attachUUID(testInputObject2)
    expect_error(logEvent(eventTitle, eventCall, input = list(testInputObject1, testInputObject2)), NA)

    # call function with uuid which does not occur in the logs
    # expect function to not return any instances
    searchUUID <- "tc17e152a-e0f2-4b9c-a70e-40ed417dc0f7"
    expect_equal(findUUID(searchUUID), NULL)
    expect_true(file.remove(logFile))
})

test_that("one instance is returned if uuid occurs in the logs once", {
    # set up log file with one instance of a uuid which will be called by the function
    # new logFile is created in tests folder
    logFile <- logFileName(setOption = TRUE)
    expect_false(file.exists(logFile))

    eventTitle <-"Graph generation"
    eventCall <- "function1(param1, param2)"

    testInputObject1 <- c("testparam1", "testparam2", "testparam3")
    attr(testInputObject1, "testBlankAttr1") <- "testBlankValue1"
    testInputObject1 <- attachUUID(testInputObject1)

    testInputObject2 <- c("testparam4")
    attr(testInputObject2, "testBlankAttr2") <- "testBlankValue2"
    testInputObject2 <- attachUUID(testInputObject2)
    expect_error(logEvent(eventTitle, eventCall, input = list(testInputObject1, testInputObject2)), NA)

    # print(readChar(logFile, file.info(logFile)$size))

    # call function with same uuid
    searchUUID <- attr(testInputObject1, "UUID")

    # expect function to return one instance
    expectedResult0 <- "Occurrence 1:\n"
    expectedResult1 <- "event\t|\tinput\t|\tattribute\t|\ttestBlankAttr1\t|\ttestBlankValue1\n"
    expectedResult2 <- paste("event\t|\tinput\t|\tattribute\t|\tUUID\t|\t", attr(testInputObject1, "UUID"), "\n", sep = "")
    expectedResult3 <- "event\t|\tinput\t|\tattribute\t|\ttestBlankAttr2\t|\ttestBlankValue2\n"
    expectedResult4 <- paste("event\t|\tinput\t|\tattribute\t|\tUUID\t|\t", attr(testInputObject2, "UUID"), "\n", sep = "")
    expectedResult5 <- "event\t|\tcall\t|\tfunction1(param1, param2)\n"
    expectedResult6 <- paste("comment\t|\ttitle\t|\tGraph generation\ncomment\t|\tdateTime\t|\t", Sys.time(), "\n\n", sep = "")

    expect_equal(findUUID(searchUUID), paste(expectedResult0, expectedResult1, expectedResult2, expectedResult3,
                                             expectedResult4, expectedResult5, expectedResult6, sep = ""))
    expect_true(file.remove(logFile))

})

test_that("multiple instances are returned if uuid occurs in the logs more than once", {
    # set up log file with 3 instances of a uuid which will be called by the function
    logFile <- logFileName(setOption = TRUE)
    expect_false(file.exists(logFile))

    # first event which only uses testInputObject1 by reading from it
    eventTitle1 <-"Graph generation"
    eventCall1 <- "function1(param1, param2)"

    testInputObject1 <- c("testparam1", "testparam2", "testparam3")
    attr(testInputObject1, "testBlankAttr1") <- "testBlankValue1"
    testInputObject1 <- attachUUID(testInputObject1)

    testInputObject2 <- c("testparam4")
    attr(testInputObject2, "testBlankAttr2") <- "testBlankValue2"
    testInputObject2 <- attachUUID(testInputObject2)
    expect_error(logEvent(eventTitle1, eventCall1, input = list(testInputObject1, testInputObject2)), NA)


    # second event which only uses testInputObject1 by reading from it
    eventTitle2 <- "More graph generation"
    eventCall2 <- "function2(param1, param2)"

    testOutputObject2 <- "testparam5"
    testOutputObject2 <- attachUUID(testOutputObject2)
    expect_error(logEvent(eventTitle2, eventCall2, input = list(testInputObject1), output = list(testOutputObject2)), NA)

    # third event which only uses testInputObject1 by reading from it
    eventTitle3 <- "Do something with graph"
    eventCall3 <- "function3(param1, param2)"

    testOutputObject3 <- "testparam6"
    testOutputObject3 <- attachUUID(testOutputObject3)
    expect_error(logEvent(eventTitle3, eventCall3, input = list(testInputObject1), output = list(testOutputObject3)), NA)

    # call function with same uuid
    searchUUID <- attr(testInputObject1, "UUID")

    # expect function to return one instance

    expectedResult0 <- "Occurrence 1:\n"
    expectedResult1 <- "event\t|\tinput\t|\tattribute\t|\ttestBlankAttr1\t|\ttestBlankValue1\n"
    expectedResult2 <- paste("event\t|\tinput\t|\tattribute\t|\tUUID\t|\t", attr(testInputObject1, "UUID"), "\n", sep = "")
    expectedResult3 <- "event\t|\tinput\t|\tattribute\t|\ttestBlankAttr2\t|\ttestBlankValue2\n"
    expectedResult4 <- paste("event\t|\tinput\t|\tattribute\t|\tUUID\t|\t", attr(testInputObject2, "UUID"), "\n", sep = "")
    expectedResult5 <- "event\t|\tcall\t|\tfunction1(param1, param2)\n"
    expectedResult6 <- paste("comment\t|\ttitle\t|\tGraph generation\ncomment\t|\tdateTime\t|\t", Sys.time(), "\n\n", sep = "")

    expectedResult7 <- "Occurrence 2:\n"
    expectedResult8 <- "event\t|\tinput\t|\tattribute\t|\ttestBlankAttr1\t|\ttestBlankValue1\n"
    expectedResult9 <- paste("event\t|\tinput\t|\tattribute\t|\tUUID\t|\t", attr(testInputObject1, "UUID"), "\n", sep = "")
    expectedResult10 <- paste("event\t|\toutput\t|\tattribute\t|\tUUID\t|\t", attr(testOutputObject2, "UUID"), "\n", sep = "")
    expectedResult11 <- "event\t|\tcall\t|\tfunction2(param1, param2)\n"
    expectedResult12 <- paste("comment\t|\ttitle\t|\tMore graph generation\ncomment\t|\tdateTime\t|\t", Sys.time(), "\n\n", sep = "")

    expectedResult13 <- "Occurrence 3:\n"
    expectedResult14 <- "event\t|\tinput\t|\tattribute\t|\ttestBlankAttr1\t|\ttestBlankValue1\n"
    expectedResult15 <- paste("event\t|\tinput\t|\tattribute\t|\tUUID\t|\t", attr(testInputObject1, "UUID"), "\n", sep = "")
    expectedResult16<- paste("event\t|\toutput\t|\tattribute\t|\tUUID\t|\t", attr(testOutputObject3, "UUID"), "\n", sep = "")
    expectedResult17 <- "event\t|\tcall\t|\tfunction3(param1, param2)\n"
    expectedResult18 <- paste("comment\t|\ttitle\t|\tDo something with graph\ncomment\t|\tdateTime\t|\t", Sys.time(), "\n\n", sep = "")

    expect_equal(findUUID(searchUUID), paste(expectedResult0, expectedResult1, expectedResult2, expectedResult3,
                                             expectedResult4, expectedResult5, expectedResult6, expectedResult7,
                                             expectedResult8, expectedResult9, expectedResult10, expectedResult11,
                                             expectedResult12, expectedResult13, expectedResult14, expectedResult15,
                                             expectedResult16, expectedResult17, expectedResult18, sep = ""))
    expect_true(file.remove(logFile))

    # call function with same uuid
    # expect function to return the 3 instances of uuid being mentioned, in chronological order
})

# tests for: getProvenance(uuid, uuidPath)

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
