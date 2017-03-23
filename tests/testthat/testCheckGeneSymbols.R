# testDummy.R
#
# Dummy tests as a template for unit testing in test-driven development
# for the rete project. If you test code that writes information to log-files,
# turn logging off for all tests except for one specific block in which
# you verify and validate your logs. Do not use your own log-file but
# rely on the global option taht is set in the setup block. We have good
# code samples for most test you will need - if you are not sure where to find
# them, ask. Do not reinvent the wheel, but write your code according to the
# principles you find in the existing files.
#
# Structure your tests according to the three principles expressed below:
#    test_that("parameter errors are correctly handled", { ... })
#    test_that("a sane input gives an expected output", { ... })
#    test_that("a corrupt input does not lead to corrupted output", { ...})


context("checkGeneSymbols functions")


# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()
# ==== END SETUP AND PREPARE ===================================================

test_that("isGeneSymbol identifies invalid gene symbols", {
    # Input contains invalid (non-HGNC) gene symbols.
    input <- c("", "!@#", "AR3bsdt!", "123")
    expected_output <- c(FALSE, FALSE, FALSE, FALSE)
    expect_equal(isGeneSymbol(input), expected_output)

})


test_that("isGeneSymbol identifies correct gene symbols, and input is case insensitive", {
    # Input contains correct HGNC gene symbols.
    input <- c("A1BG", "a1bg", "A1bG", "a1Bg")
    expected_output <- c(TRUE, TRUE, TRUE, TRUE)
    expect_equal(isGeneSymbol(input), expected_output)
})


test_that("a corrupt input does not lead to corrupted output", {
    # Input is corrupted (NULL, NA, vector of length zero, or no input at all).
    expect_error(isGeneSymbol())
    expect_error(isGeneSymbol(c()))
    expect_error(isGeneSymbol(NULL))
    expect_error(isGeneSymbol(NA))
})



# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
