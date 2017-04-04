#Pseudocode may look like actual code...please ignore...
get_pdb_ids <- function(target_seq, start_pos = 1 , end_pos = nchar(target_seq), convert = TRUE, silent = TRUE){
  
  #get packages
  # install.packages("bio3d")
  # install.packages("RCurl")
  # install.packages("XML")
  # install.packages("seqinr")
  # library(bio3d)
  # library(RCurl)
  # library(XML)
  # library(seqinr)
  
  
  
  #convert sequence into AA 
  if(convert){
    #translate needs a vector of single chars, not a string
    AA_seq <- paste(translate(s2c(target_seq)), collapse='')
    end_pos <- nchar(AA_seq)
  }else{ AA_seq <- target_seq }
 
  blast_output <- blast.pdb(AA_seq)
  
  blast_output <- blast_output$hit.tbl
  blast_output <- blast_output[blast_output$identity > 40, ]
  
  if(nrow(blast_output) == 0){
    return(NULL)
  }
  
  blast_output$tax <- NA 
  blast_output$info <- NA
  
  final_output <- data.frame()
  final_output_range <- c()
  final_output_range_counter <- 1
  
  for(entry in 1:nrow(blast_output)){
    rcsb_id <- gsub("_", ".", blast_output$pdb.id[entry], ignore.case = TRUE)
    lookup <- xmlParse(getURL(paste("http://www.rcsb.org/pdb/rest/describeMol?structureId=", rcsb_id, sep = "")))
    
    #get taxonomy 
    
    tax <- getNodeSet(lookup, '//Taxonomy/@name')[[1]][[1]]
    info <- getNodeSet(lookup, '//polymerDescription/@description')[[1]][[1]]
    
    blast_output$tax[entry] <- ifelse(length(tax) > 0, tax, NA)
    blast_output$info[entry] <- ifelse(length(info) > 0, yes = info, no = NA)
    
    if(!intersecting_positions(final_output_range, blast_output$q.start[entry], blast_output$q.end[entry])){
      if((!is.na(blast_output$tax[entry]) && blast_output$tax[entry] == "Homo sapiens") || blast_output$identity[entry] == 100){
        final_output <- rbind(final_output, blast_output[entry,])
        final_output_range[[final_output_range_counter]] <- c(blast_output$q.start[entry]:blast_output$q.end[entry])
        final_output_range_counter <- final_output_range_counter + 1
      }
    }
  }
  
  if(nrow(final_output) == 0){ return(NULL) }
  
  final_output$q.start <- final_output$q.start + start_pos - 1
  final_output$q.end <- final_output$q.end + start_pos - 1
  
  for(x in 1:length(final_output_range)){
    final_output_range[[x]] <- final_output_range[[x]] + start_pos - 1
  }
  
  min_side <- min(final_output_range[[1]])
  max_side <- max(final_output_range[[length(final_output_range)]])
  
  if(end_pos - max_side > 50){
    message(sprintf("running substr from pos %d to %d", max_side, end_pos))
    max_frame <- get_pdb_ids(substr(AA_seq, max_side - start_pos, end_pos - start_pos), 
                             start_pos = max_side, end_pos = end_pos, convert = FALSE, silent = silent)
    final_output <- rbind(final_output, max_frame)
  }
  
  if(min_side > 50){
    message(sprintf("running substr from pos %d to %d", start_pos, min_side))
    min_frame <- get_pdb_ids(substr(AA_seq, 1, min_side - start_pos), 
                             start_pos = start_pos, end_pos = min_side, convert = FALSE, silent = silent)
    final_output <- rbind(final_output, min_frame)
  }
  
  for(y in 1:(length(final_output_range) - 1)){
    left_side <- max(final_output_range[[y]])
    right_side <- min(final_output_range[[y+1]])
    
    if(right_side - left_side > 50){
      message(sprintf("running substr from pos %d to %d", left_side, right_side))
      inter_frame <- get_pdb_ids(substr(AA_seq, left_side - start_pos, right_side - start_pos), 
                               start_pos = left_side, end_pos = right_side, convert = FALSE, silent = silent)
      final_output <- rbind(final_output, inter_frame)
    }
  }
  
  return(final_output)
}

intersecting_positions <- function(r, st, ed){
  for(x in 1:length(r)){
    if(length(intersect(r[[x]], c(st:ed))) > 0){
      return(TRUE)
    }
  }
  return(FALSE)
}