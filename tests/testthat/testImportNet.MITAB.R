#Check fName input is correct
expect_error(fName = NULL || fName = c("a", "b", "c", "d") || "nosuchtext.txt" ||, raise error)
#Check drop unmapped
expect_error(dropUnmapped = NULL, raise error)
#Check silent
expect_error(silent = NULL, raise error)
#chek write log
expect_error(writeLog = NULL, rasie error)

test_that(importNet.MITAB creates gG from MITAB data){
  ##fName containing only 5 valid vertices... I dont know how this will be done?
  expect_equal(igraph::vcount(gG), 5)
  expect_equal(igraph::ecount(gG), 4)
  
  #Check correct meta data
  check gG type
  check gG version
  check UUID
  
  #Check log file
  Check if log file created correctly
  
  #check calls
  check call
  check fname
  check val
  check cutofftype
  check taxID
  check dropUnmapped
  check silent
  check writelog
  check event note
  check edgges selectM
  check gG vertices and edges
  check gG
  check event output
  check class
  check type
  
  #check leaks
  check no loak output if silent
  
  
}
