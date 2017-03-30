importNet.STRING <- function(fName,
                             cutoffType = "xN",
                             val,
                             experimentType = getOptions("rete.defaultPPI"),
                             taxID = "9606",
                             silent = FALSE,
                             writeLog = FALSE) {
  
  #Parameters
  1) set cutoffTypes
  2) set cutofftype default values
  3) use val given or default if not given
  
  #Validatations
  4) Validate inputs

  5) Validate fName input #Needed because columns "a" "b" and "weight" are not same format for all
        if iRefWeb calculate score method 1 and gsub method 1
        if inBioMap calculate score method 2 and gsub method 2
        if IntAct calculate score method 3 and gsub method 3

  6) Read file
  7) Name column 1-"a" 2-"b" and 15-"weight"
  8) Calculate score and save score into numerical
  9) Remove all rows that dont contain taxid 9606
  10) subset dataframe by methods removeByMethods()
  11) Remove all rows that dont pass the cutoff
  12) fastMap()
  13) dropunmapped()
  14) simplify()
  15) .df2gG()
  16) Writelog()
  17) Return(gG)
  
  
  
}