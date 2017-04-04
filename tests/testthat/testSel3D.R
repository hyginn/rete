# testSel3D.R
#
#
# Structure your tests according to the three principles expressed below:
#    test_that("parameter errors are correctly handled", { ... })
#    test_that("a sane input gives an expected output", { ... })
#    test_that("a corrupt input does not lead to corrupted output", { ...})


# ==== BEGIN SETUP AND PREPARE =================================================
#OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
#logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
#NL <- .PlatformLineBreak()
# ==== END SETUP AND PREPARE ===================================================




#test_that("parameter errors are correctly handled", {
  # Try missing input sequence
  # Try extraneous parameters
  # Note: we're not to attempt to avoid all kinds of abuse from users
  #   but does that include messing with the preset parameters 
  #   which only exist for the sake of recursion?
  #expect_equal()

#})


#test_that("a sane input gives an expected output", {
  # Note: All of the expected results need to be verified by BLAST to PDB
  # Try single codon
  # Try an input which has no results (verify with BLAST to PDB)
  # Try small input with a small result
  # Try a large input with many resulting matches (verify with BLAST to PDB)
#})


#test_that("a corrupt input does not lead to corrupted output", {
  # Try illegal characters in input sequence
  # Try illegal length of input sequence
  # Try empty string as input sequence
  # Try unexpected type of input sequence
  # These should all exit immediately without executing any of the actual algo.
#})



# ==== BEGIN TEARDOWN AND RESTORE ==============================================
#logName <- unlist(getOption("rete.logfile"))
#if (file.exists(logName)) { file.remove(logName)}
#options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
