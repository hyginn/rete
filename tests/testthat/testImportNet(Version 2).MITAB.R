OLOG <- as.character(getOption("rete.logfile"))
logFileName(fPath = tempdir(), setOption = TRUE)
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) {file.remove(logName)}

# ==== END SETUPS AND PREPARE ===================

# ==== importNet.MITAB ==========================

test_that("importNet.MITAB() rejects erroneus parameters", {
  fN <- "core.psimitab"
  
  # === fName
  
  expect_error(
    gG <- importNet.MITAB(fName = NULL,
                          cutoffType = "xN",
                          silent = TRUE),
    '.checkArgs> "fName" mode error')

  expect_error(
    gG <- importNet.MITAB(fName = c("file_A", "file_B"),
                          cutoffType = "xN",
                          silent = TRUE),
    '.checkArgs> "fName" length error')
  
  expect_error(
    gG <- importNet.MITAB(fName = "noSuch.txt", 
                          cutoffType = "xN",
                          silent = TRUE),
    '.checkArgs> "fName" error: file "noSuch.txt" does not exist.')
  
  # ==== dropUnmapped
  
  expect_error(
    gG <- importNet.MITAB(fName = fN,
                          cutffType = "xN",
                          dropUnmapped = NULL,
                          silent = TRUE),
    '.checkArgs> "dropUnmapped" mode error')
  
  expect_error(
    gG <- importNet.MITAB(fName = fN,
                          cutoffType = "xN",
                          silent = NULL),
    '.checkArgs> "silent" mode error')
  
  expect_error(
    gG <- importNet.MITAB(fNAme = fN,
                          cutoffType = "xN",
                          silent = TRUE,
                          writeLog = NULL),
    '.checkArgs> "writeLog" mode error')
  
  expect_false(file.exists(logName))

}




importNet.STRING <- function(fName,
                             cutoffType = "xN",
                             val,
                             experimentType = getOptions("rete.defaultPPI"),
                             taxID = "9606",
                             silent = FALSE,
                             writeLog = FALSE) {

})

test_that("importNet.MITAB() correctly calculates scores", {
  fN <- "core.psimitab"
  g
})



test_that("importNet.MITAB() produces gG from MITAB data", {
  unlink(logName) # Deleting previous log file?
  fN <- "core.psimitab"
  gG <- importNet.MITAB(fName = fN,
                        cuttoffType = "xN",
                        dropUnmapped = FALSE,
                        silent = TRUE,
                        weriteLog = TRUE)
  expect_equal(igraph::vcount(gG), 5)
  expect_equal(igraph::ecount(gG), 4)
  
  expect_equal(attr(gG, "type"), "gG")
  expect_equal(attr(gG, "version"), "1.0")
  expect_equal(.checkArgs(attr(gG, "UUID"), "UUID"), character())
  
  thisLog <- readlines(logName)
  expect_true(grepl("event \\| title \\| importNet.MITAB",     thisLog[1]))
  expect_true(grepl("event \\| time",                          thisLog[2]))
  
  expect_true(grepl("^event \\| call \\| importNet.MITAB \\(", thisLog[3]))
  expect_true(grepl(paste0("fName = \"", fN, "\""),            thisLog[3]))
  expect_true(grepl("cutoffType = \"xN\"",                     thisLog[3]))
  expect_true(grepl("val = 10000",                             thisLog[3]))
  expect_true(grepl("taxID = \"9606\"",                        thisLog[3]))
  expect_true(grepl("dropUnmapped = FALSE",                    thisLog[3]))
  expect_true(grepl("silent = TRUE",                           thisLog[3]))
  expect_true(grepl("writeLog = TRUE",                         thisLog[3]))
  expect_true(grepl(")$",                                      thisLog[3]))
  
  expect_true(grepl("event | note | read 4 edges from file.",  thisLog[4]))
  expect_true(grepl("Selected 4 edges via cutooff.",           thisLog[5]))
  expect_true(grepl("gG object has 5 vertices and 4 edges.",   thisLog[6]))
  expect_true(grepl("\"gG\"",                                  thisLog[7]))
  expect_true(all(grepl("event \| output \\| attribute \\|", thisLog[8:11])))
  expect_true(grepl("class | \"igraph\"",                      thisLog[8]))
  expect_true(grepl("type | \"gG\"",                           thisLog[9]))
  expect_true(grepl("version | \"1.0\"",                       thisLog[10]))
  expect_true(grepl("UUID",                                    thisLog[11]))
  expect_true(grepl("event | end",                             thisLog[12]))
  expect_true(grepl("^$",                                      thisLog[13]))
  })

  test_that("importNet.MITAB() does not leak output if silent = TRUE", {
    fN <- "core.psimitab"
    testF <- tempfile()
    capture.output(gG <- importNet.MITAB(fName = fN,
                                         cutoffType = "xN",
                                         dropUnmapped = FALSE,
                                         silent = TRUE,
                                         writeLog = FALSE),
                   file = testF)
    expect_equal(length(readLines(testF)), 0)
    expect_true(file.remove(testF))
  })
  
  
# ==== LOOK HERE ============================================================================  
  test_that("importNet.MITAB() follows the designated cutoff", {
    fN <- system.file("extdata", "core_snippet.psimitab", package="rete")
    
    val <- 5 # Top 5 highest scores
    gG <- importNet.MITAB(fName = fN,
                          cutoffType = "xN",
                          dropUnmapped = TRUE,
                          silent = TRUE,
                          writeLog = FALSE)
    
    expect_true(gorder(gG), 5) # Check quantity of vertices 
    expect_true(E(gG)$weights[1], 0.190986) # Check weight of first edge
    
    
    val <- 0.45 # Edges fall into 0.45 quartile
    gG <- importNet.MITAB(fName = fN,
                          cutoffType = "xQ",
                          val,
                          dropUnmapped = TRUE,
                          silent = TRUE,
                          writeLog = FALSE)
    expect_true(gorder(gG), 5) # Check quantity of vertices 
    
    val <- 0.25 # Edges with weight >= 0.25
    gG <- importNet.MITAB(fName = fN,
                          cutoffType = "xS",
                          val,
                          dropUnmapped = TRUE,
                          silent = TRUE,
                          writeLog = FALSE)
    expect_true(gorder(gG), 5) # Check quantity of vertices 
  })
  
# ==== HERE ================================================================================
  
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
  
  
  
  
  
  
  
  
  
  
  
  