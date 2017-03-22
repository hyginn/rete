# testCombine.R

context("test for Combine SNV and CNA files")

# ==== BEGIN SETUP AND PREPARE =================================================
# set up tempdir and tempfiles for output files and to be deleted after tests
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
# ==== END SETUP AND PREPARE ===================================================

test_that("parameter errors are correctly handled", {
    # try no parameter input
    expect_error(combineSNV_CNA(writeLog=FALSE))
    # try length-zero or NuLL input vector
    input_vector <- character()
    expect_error(combineSNV_CNA(input_vector, writeLog=FALSE))
    # try wrong file format (e.g. text file instead of rds dataframe) in input vector
    input_vector[1] <- "CNA.txt"
    expect_error(combineSNV_CNA(input_vector, writeLog=FALSE))
    # try invalid output fgX file name (e.g. wrong extension)
    input_vector[1] <- "CNA.rds"
    expect_error(combineSNV_CNA(input_vector, "out.txt", writeLog=FALSE))
})


test_that("a sane input gives an expected output", {
    # provide one valid SNV and one valid CNA file with small number of genes as input
    # TODO where can I find the sample data for testing? inst/extdata folder has maf & txt formats
    input_vector = c("CNA.rds", "SNV.rds")

    # test if a log file is saved if log is True
    combineSNV_CNA(input_vector, writeLog=TRUE)
    expect_true(files.exists(logName))
    # test if output file fgX is saved under the path
    combineSNV_CNA(input_vector, writeLog=FALSE)
    expect_true(files.exists("gx.rds"))
    # test if specified output filename is saved under the path
    combineSNV_CNA(input_vector, fgx = "test.rds", writeLog=FALSE)
    expect_true(files.exists("gx.rds"))

})


test_that("combine received corrupted input data file", {
    # try SNV file with missing data (e.g. only gene symbol, no mutation class), expect error and stop in log and console
    # TODO edit a CNA file to which some headers have empty/null data called corrupted_CNA.rds
    input_vector = c("corrupted_CNA.rds", "SNV.rds")
    expect_error(combineSNV_CNA(input_vector, writeLog=FALSE))
    # try CNA file with missing copynumber, expect error and stop in log and console
    # TODO edit a SNV file to which gene symbols correponds to null variant class & position called corrupted_SNV.rds
    input_vector = c("CNA.rds", "corrupted_SNV.rds")
    expect_error(combineSNV_CNA(input_vector, writeLog=FALSE))
})

# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================
# [END]
