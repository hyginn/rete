library(stringr)

importNet.STRING <- function(fName,
                             cutoffType = "xN",
                             val,
                             experimentType = getOptions("rete.defaultPPI"),
                             taxID = "9606",
                             silent = FALSE,
                             noLog = FALSE) {
  
  # ==== PARAMETERS =========================================================
  cutoffTytpes <- c("xS", "xQ", "xN")
  defaultValues <- c(xS = 950, xQ = 0.9, xN = 10000)
  
  if (missing(val)) {
    val <- defaultValues[cutoffType]
  }
  
  # ==== VALIDATIONS =======================================================
  
  cR <- character()
  cR <- c(cR, .checkArgs(fName,        like = "FILE_E",   checkSize = TRUE))
  cR <- c(cR, .checkArgs(cutoffType,   like = "a",        checkSize = TRUE))
  cR <- c(cR, .checkArgs(val,          like = numeric(1), checkSize = TRUE))
  cR <- c(cR, .checkArgs(taxID,        like = "a",        checkSize = TRUE))
  cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
  cR <- c(cR, .checkArgs(noLog,        like = logical(1), checkSize = TRUE))
  
  if(length(cR) > 0) {
    stop(cR)
  }
  
  # Validate the cutoffType parameter
  if (!cutoffType %in% cutoffTypes) {
    stop(paste0("Parameter error:", "\n",
                "Valid cutoff types: ",
                "\"", paste(cutoffTypes, collapse = "\" | \""), "\"", "\n",
                "Found: \"", cutoffType, "\"", "\n"))
  }
  
  # Read data and create header
  
  readfile <- read.delim("core.psimitab", sep = "\t", header = FALSE)
  colnames(readfile) <- c("a", "b", "c", "d","e","r","t","y","u","i","o","p","[","s","weight","f" )
  scores <- readfile$weight
  
  # Create single score without | operator
  i = 1
  weightlist <- NULL
  while(i <= length(scores)){
    variable <- unlist(strsplit(as.character(scores[i]), "\\|"))
    if(variable[1] != "-" && variable[2] != "-"){
      value <- as.numeric(variable[1]) * as.numeric(variable[2])
      weightlist[i] <- value
    } else weightlist[i] <- 1
    i = i + 1 
  }
  # Simplify data into only the needed coloumns 
  readfile$weight <- weightlist
  a <- which(colnames(readfile) == "a")
  b <- which(colnames(readfile) == "b")
  c <- which(colnames(readfile) == "weight")
  netStr <- readfile[c(a,b,c)]
  
  if (cutoffType == "xS") {
    # select all rows with weight >= val
    sel <- netStr$weight >= val
    nSel <- sum(sel)
  } else if (cutoffType == "xQ"){
    x <- stats::quantile(netDF$weight, probs = val)
    sel <- netStr$weight <= val
    nSel <- sum(sel) 
  } else if (cutoffType == "xN") {
    # select the val highest scores, or all, whichever is fewer
    x <- order(netStr$weight, decreasing = TRUE)
    sel <- x[1:min(val, length(x))]
    nSel <- min(val, length(x))
  }
  readfile <- readfile[sel, ] # Keep only selected genes
  
  fastMapENSP <- readRDS(system.file("extdata", "fastMapENSP.rds", package = "rete"))
  
  #Remove prefix from PSIMITAB format IDs
  netStr$a <- gsub(".*:", "", netStr$a)
  netStr$b <- gsub(".*:", "", netStr$b)
  
  # map ID to gene names
  mapA <- fastMap(netStr$a, fastMapENSP, type = "ENSP", dev = TRUE)
  mapB <- fastMap(netStr$b, fastMapENSP, type = "ENSP", dev = TRUE)
  
  if (dropUnmapped) {
    netStr$a <- mapA
    netStr$b <- mapB
    # select edges where one or the other gene symbol is NA
    sel <- is.na(netStr$a) | is.na(netStr$b)
    nDrop <- sum(sel)
    netDF <- netStr[-sel, ]  # Drop unmapped edges
  } else {
    # replace mapped IDs with their gene symbol
    netStr$a[! is.na(mapA)] <- mapA[! is.na(mapA)]
    netStr$b[! is.na(mapB)] <- mapB[! is.na(mapB)]
  }
  
  gG <- .df2gG(netStr)
  
}

