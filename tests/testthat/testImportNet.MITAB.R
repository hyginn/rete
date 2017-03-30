#bs> rememeber to add setup/teardown blocks to the file

#Check fName input is correct
expect_error(fName = NULL || fName = c("a", "b", "c", "d") || "nosuchtext.txt" ||, raise error)
#bs> don't need to check 4 if you want to check > 1
#Check drop unmapped
expect_error(dropUnmapped = NULL, raise error)
#Check silent
expect_error(silent = NULL, raise error)
#chek write log
expect_error(writeLog = NULL, rasie error)
#bs> check that after all ht efailed parameter tests no log-file was created
#bs> remember to confirm the actual error message, not just that an error was thrown





test_that(importNet.MITAB creates gG from MITAB data){
  ##fName containing only 5 valid vertices... I dont know how this will be done?
#bs> Can you grep out a few records relating to our test-dataset?
  
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

#bs> check that calculated scores are correct
#bs> check that selections by value / number quartile ...  work
#bs> check that ID translation works
    
#bs> check that unmapped IDs are correctly dropped if requested
#bs> check that symmetric edges are correctly added if requested
#bs> check that edges that don't match requested methods are properly removed

#bs> check that number of unmapped IDs was correctly reported
#bs> check that edges added as symmetric edges are properly reported
#bs> check that edges removed via method selection are properly reported

#bs> check that resulting graph is weighted, directed and simple
    
    
    
  
}
