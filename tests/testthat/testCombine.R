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
                 'Input vectors are empty')

    # try wrong file format (e.g. text file instead of rds dataframe) in SNV input vector
    input_vector = c("SNV.txt", "CNA.rds")
    expect_error(combineSNV_CNA(input_vector[1], input_vector[2], writeLog=FALSE),
                 'SNV file has to be .rds format.')

    # try wrong file format (e.g. text file instead of rds dataframe) in CNA input vector
    input_vector = c("SNV.rds", "CNA.txt")
    expect_error(combineSNV_CNA(input_vector[1], input_vector[2], writeLog=FALSE),
                 'CNA file has to be .rds format.')

    # sane input vector
    input_vector = c("SNV.rds", "CNA.rds")
    # try non existing input file for SNV input vector
    expect_error(combineSNV_CNA(c("nonExistent.rds"), input_vector[2], writeLog=FALSE),
                 '.checkArgs> "fSNV" error: file "nonExistent.rds" does not exist')

    # try non existing input file for CNA input vector
    expect_error(combineSNV_CNA(input_vector[1], c("nonExistent.rds"), writeLog=FALSE),
                 '.checkArgs> "fCNA" error: file "nonExistent.rds" does not exist')

    # try invalid output fgX file name (e.g. wrong extension)
    expect_error(combineSNV_CNA(input_vector[1], input_vector[2], "out.txt", writeLog=FALSE),
                 'Output file has to be .rds extension!')

    # try null output fgX file name
    expect_error(combineSNV_CNA(input_vector[1], input_vector[2], NULL, writeLog=FALSE),
                 'argument is of length zero')

    # try null silent option
    expect_error(combineSNV_CNA(input_vector[1], input_vector[2], silent=NULL, writeLog=FALSE),
                 '.checkArgs> "silent" mode error')

    # try null writeLog option
    expect_error(combineSNV_CNA(input_vector[1], input_vector[2], writeLog=NULL),
                 '.checkArgs> "writeLog" mode error')

    # check if log file has been written
    expect_false(file.exists(logName))
})


test_that("a sane input gives an expected output", {
    # provide one valid SNV and one valid CNA file with small number of genes as input
    # TODO where can I find the sample data for testing? inst/extdata folder has maf & txt formats
    input_vector = c("CNA.rds", "SNV.rds")
    testF <- tempfile(fileext = ".rds")

    # test if a log file is saved if log is True
    capture.output(combineSNV_CNA(input_vector[1], input_vector[2], fgX=testF, silent=TRUE, writeLog=TRUE), file=testF)
    expect_true(files.exists(logName))

    # correct log file
    thisLog <- readLines(logName)
    expect_true(grepl("event \\| title  \\| combineSNV_CNA",       thisLog[1]))
    expect_true(grepl("event \\| time",                            thisLog[2]))

    # call
    expect_true(grepl("^event \\| call   \\| combineSNV_CNA\\(",   thisLog[3]))
    expect_true(grepl(paste0("fNameSNV = \"", input_vector[1], "\""),    thisLog[3]))
    expect_true(grepl(paste0("fNameCNA = \"", input_vector[2], "\""),    thisLog[3]))
    expect_true(grepl("fgX = gX.rds",                              thisLog[3]))
    expect_true(grepl("silent = TRUE",                             thisLog[3]))
    expect_true(grepl("writeLog = TRUE",                           thisLog[3]))
    expect_true(grepl(")$",                                        thisLog[3]))

    expect_true(grepl("event | note   | Combine 10 samples in total.",  thisLog[4]))
    expect_true(grepl("event | end",                               thisLog[5]))
    expect_true(grepl("^$",                                        thisLog[6]))


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
