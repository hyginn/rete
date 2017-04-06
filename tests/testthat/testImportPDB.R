#testImportPDB.R

context("<Dummy tests>")


# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()
# ==== END SETUP AND PREPARE ===================================================

test_that("parameter errors are correctly handled", {
    # Try each parameter missing in turn
    # Try parameters out of expected bounds
    # Try NULL and length-zero parameters: be sure none lead to an erroneous
    #    one-time execution of any loop
    expect_error(importPDB())
    expect_error(importPDB(TRUE))
    expect_error(importPDB("2"))
    expect_error(importPDB(NULL))
    expect_error(importPDB(numeric()))
})


test_that("a sane input gives an expected output", {
    # Provide small input data, provide small output data.
    #   (Note that the output needs to be independently verifiably correct.
    #    Don't run an erroneous function, and then test against the erroneous
    #    output you received when you ran your function for the first time.)
    # Try this with important variations of parameters (reasonable, not
    #    combinatorially exhaustive).
    # Cover your code.

    #TODO...
})


test_that("a corrupt input does not lead to corrupted output", {
    # Try: - spurious characters in numeric columns ("N/A", ...).
    #      - extra tabs at line-end
    #      - comments before headers
    #      - truncated file (incomplete last line)
    # Again: make absolutely sure be you never have an erroneous
    #        one-time execution of a loop
    expect_error(importPDB(NA))
})



# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
