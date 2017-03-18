# testImportSNV.TCGA.R

context("Tools for importing MAF files into a rSNV file")

test_that("parameter errors are correctly handled", {
    # fMAF is not provided
    # fMAF is NULL
    # fMAF input vector is empty
    # fMAF input vector contains non-string values
    # file in fMAF could not be opened/does not exist
    # if rSNV parameter is provided, rSNV file does not exist
    # if rSNV parameter is provided, but is NULL
    # if rSNV parameter is provided, but is non-string value
    # if rSNV is not provided, new file could not be created (in case of no  write permissons)
    # if na.rm is not given logical input
    # if silent is not given logical input
    # write.log is not given logical input
})

test_that("check various forms of input parameters work correctly", {
    # only fMAF is provided
    # fMAF, rSNV provided
    # fMAF, na.rm
    # fMAF, silent
    # fMAF, writeLog
})


test_that("a sane input gives an expected output", {
    # make small input
    # make small output for that input manually
    # run the module on input
    # check if manual output = module output
})


test_that("a corrupt input does not lead to corrupted output", {
    # if fMAF file open does not follow the 2.4 or agreed upon version either STOP or DROP file
    # if fMAF file open does not have the required columns at the
    #           right column number STOP or DROP file
    # if rSNV is provided, check if number of columns in rSNV is worng
    # if rSNV is provided, check when header labels are wrong
})

test_that("na.rm works as intended", {
    # if fMAF file has N/A values in rows. if na.rm == TRUE. check if the row
    #           is dropped and not included in rSNV
    # if fMAF file has N/A values in rows. if na.rm == FALSE check if the row
    #           is not dropped and included in rSNV

})

test_that("silent and write works as intended", {
    # if silent=TRUE, check if output to console is supressed
    # if silent=FALSE, check if output to console is done correctly
    # if writeLog=TRUE, check if log file is created
    # if writeLog=TRUE, check if logs are written to it correctly
    # if writeLog=FALSE, check that no log file is created
    # if writeLog=FALSE, cmake sure no log files are written
})

# [END]
