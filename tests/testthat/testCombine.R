# testCombine.R

context("test for Combine SNV and CNA files")

# ==== BEGIN SETUP AND PREPARE =================================================
# set up tempdir and tempfiles for output files and to be deleted after tests
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
# ==== END SETUP AND PREPARE ===================================================

test_that("parameter errors are correctly handled", {
    # try no parameter input
    expect_error(combineSNV_CNA(writeLog=FALSE),
                 '.checkArgs> "fName" mode error')

    # try length-zero or NuLL input vector
    input_vector <- character()
    expect_error(combineSNV_CNA(input_vector, writeLog=FALSE),
                 '.checkArgs> "fName" length error')

    # try wrong file format (e.g. text file instead of rds dataframe) in input vector
    input_vector[1] <- "CNA.txt"
    expect_error(combineSNV_CNA(input_vector, writeLog=FALSE),
                 '.checkArgs> "fName" mode error')

    # try non existing input file
    input_vector[1] <- "nonExistent.rds"
    expect_error(combineSNV_CNA(input_vector, writeLog=FALSE),
                 '.checkArgs> "fName" error: file "NonExistent.rds" does not exist')

    # try invalid output fgX file name (e.g. wrong extension)
    input_vector[1] <- "CNA.rds"
    expect_error(combineSNV_CNA(input_vector, "out.txt", writeLog=FALSE),
                 'invalid output file extension')

    # try null output fgX file name
    expect_error(combineSNV_CNA(input_vector, NULL, writeLog=FALSE),
                 '.checkArgs> "fgx" mode error')

    # try null silent option
    expect_error(combineSNV_CNA(input_vector, silent=NULL, writeLog=FALSE),
                 '.checkArgs> "silent" mode error')

    # try null writeLog option
    expect_error(combineSNV_CNA(input_vector, writeLog=NULL),
                 '.checkArgs> "writeLog" mode error')

    # check if log file has been written
    expect_false(file.exists(logName))
})


test_that("a sane input gives an expected output", {
    # provide one valid SNV and one valid CNA file with small number of genes as input
    # TODO where can I find the sample data for testing? inst/extdata folder has maf & txt formats
    input_vector = c("CNA.rds", "SNV.rds")
    # testdir <- tempdir()
    testF <- tempfile(fileext = ".rds")
    # test_path = file.path(testdir, testF)
    print(test_path)

    # test if a log file is saved if log is True
    capture.output(combineSNV_CNA(input_vector, fgX=testF, silent=TRUE, writeLog=TRUE), file=testF)
    expect_true(files.exists(logName))

    # correct log file
    thisLog <- readLines(logName)
    expect_true(grepl("event \\| title  \\| combineSNV_CNA",       thisLog[1]))
    expect_true(grepl("event \\| time",                            thisLog[2]))

    # call
    expect_true(grepl("^event \\| call   \\| combineSNV_CNA\\(",   thisLog[3]))
    expect_true(grepl(paste0("fname = \"", input_vector, "\""),    thisLog[3]))
    expect_true(grepl("fgX = gX.rds",                              thisLog[3]))
    expect_true(grepl("silent = TRUE",                             thisLog[3]))
    expect_true(grepl("writeLog = TRUE",                           thisLog[3]))
    expect_true(grepl(")$",                                        thisLog[3]))

    expect_true(grepl("event | end",                               thisLog[4]))
    expect_true(grepl("^$",                                        thisLog[5]))

    # test if output file fgX is saved under the path
    expect_true(files.exists(testF))
    # remove tempfile
    expect_true(file.remove(testF))

})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================
# [END]
