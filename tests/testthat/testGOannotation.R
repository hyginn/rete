# testGOannotation.R

context("test for building DAG of Gene ontology with Gene ontology annotations")

# ==== BEGIN SETUP AND PREPARE =================================================
# set up tempdir and tempfiles for output files and to be deleted after tests
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
# ==== END SETUP AND PREPARE ===================================================

test_that("parameter errors are correctly handled", {
    # try no parameter input
    # try non-existing input files for both .obo and .gaf files
    # try invalid output files for both fGDAG and fGg
    # try null output
    # try null silent
    # try null writeLog
})

test_that("testGOannotation does not proceed if roots are not ontologies") {
    # try catching error with data set having the GO term for Molecular function has parents
    # try catching error with data set having the GO term for Biological process has parents
    # try catching error with data set having the GO term for Cellular component has parents
}

test_that("a sane input gives an expected output", {
    # try input with ideal dataset described on task page
    # test correct logfile
    # test if two .RDS files are saved
    # test if is_obsolete terms are removed
    # test if the minimum distance for a leaf is correct
    # test if the number of genes annotated to GO term is correct
    # test the number of GO terms left in the DAG is as expected
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================
# [END]
