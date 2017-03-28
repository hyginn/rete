# testCombine.R
#

context("Combine SNV and CNA files")


test_that("parameter errors are correctly handled", {
    # try length-zero or NuLL input vector
    # try wrong file format (e.g. text file instead of rds dataframe) in input vector
    # try files not belong to either SNV/CNA in input vector
    # try invalid output fgX file name (e.g. already existed, wrong extension)

})


test_that("a sane input gives an expected output", {
    # provide one valid SNV and one valid CNA file with small number of genes as input
    # test if returned number is correct
    # test if a log file is saved if log is True
    # test if process are sent to console is silent is True
    # test if output file fgX is saved under the path

})


test_that("combine received corrupted input data file", {
    # try both empty SNV and CNA datafile, expect empty fgX dataframe
    # try SNV file with missing data (e.g. only gene symbol, no mutation class), expect error and stop in log and console
    # try CNA file with missing copynumber, expect error and stop in log and console

})


# [END]
