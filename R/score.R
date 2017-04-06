# score.R

#' Generate the gH hash table.
#'
#' \code{score} will evaluation mutation data for genes and assign a score.
#' The score is a likelihood that the observed mutation rate exceeds the
#' background mutation rate, assuming a normal distribution. Background mutation
#' is calculated by creating a codon to amino acid table, compiling a codon usage
#' table using the GeneData resource's gene sequences, and then simulating background
#' mutation per gene with 5 simulation runs in 10 individuals. The function returns
#' a gH hash table with gene symbols as keys and values are likelihoods (scores).
#'
#'
#' @param gX A data frame containing MUT and CNA data for each gene.
#' @param GeneData A list of gene annotations.
#' @param outputName The path to save the RDS.
#' @param silent Controls whether output to console should be suppressed. FALSE
#'   by default.
#' @param writeLog Controls whether writing the result to the global logfile is
#'   enabled. TRUE by default.
#' @return a gH hash of gene symbol and associated likelihood ratio (score).
#'
#' @family score functions
#'
#' @examples
#' \dontrun{
#' score(gX, GeneData, outputName = "gH.rds", silent = FALSE, writeLog = TRUE)
#' }
#' @export

score <- function(gX, GeneData, outputName = "gH.rds", silent = FALSE, writeLog = TRUE) {


    # COMPILE CODON TO AMINO ACID TABLE

    codons <- c("AAA", "AAC", "AAT", "AAG", "ACA", "ACC", "ACT", "ACG", "ATA", "ATC", "ATT", "ATG",
                "AGA", "AGC", "AGT", "AGG", "CAA", "CAC", "CAT", "CAG", "CCA", "CCC", "CCT", "CCG",
                "CTA", "CTC", "CTT", "CTG", "CGA", "CGC", "CGT", "CGG", "TAA", "TAC", "TAT", "TAG",
                "TCA", "TCC", "TCT", "TCG", "TTA", "TTC", "TTT", "TTG", "TGA", "TGC", "TGT", "TGG",
                "GAA", "GAC", "GAT", "GAG", "GCA", "GCC", "GCT", "GCG", "GTA", "GTC", "GTT", "GTG",
                "GGA", "GGC", "GGT", "GGG")

    aa_table <- new.env()
    for(codon in codons){
        if(codon %in% c("GCT", "GCA", "GCC", "GCG")){
            aa_table[[codon]] <- "A"
        }  else if(codon %in% c("TGC", "TGT")){
            aa_table[[codon]] <- "C"
        } else if(codon %in% c("GAC", "GAT")){
            aa_table[[codon]] <- "D"
        } else if(codon %in% c("GAA", "GAG")){
            aa_table[[codon]] <- "E"
        } else if(codon %in% c("TTC", "TTT")){
            aa_table[[codon]] <- "F"
        } else if(codon %in% c("GGA", "GGC", "GGG", "GGT")){
            aa_table[[codon]] <- "G"
        } else if(codon %in% c("CAC", "CAT")){
            aa_table[[codon]] <- "H"
        } else if(codon %in% c("ATA", "ATC", "ATT")){
            aa_table[[codon]] <- "I"
        } else if(codon %in% c("AAA", "AAG")){
            aa_table[[codon]] <- "K"
        } else if(codon %in% c("TTA", "TTG", "CTA", "CTC", "CTG", "CTT")){
            aa_table[[codon]] <- "L"
        } else if(codon %in% c("ATG")){
            aa_table[[codon]] <- "M"
        } else if(codon %in% c("AAC", "AAT")){
            aa_table[[codon]] <- "N"
        } else if(codon %in% c("CCA", "CCC", "CCG", "CCT")){
            aa_table[[codon]] <- "P"
        } else if(codon %in% c("CAA", "CAG")){
            aa_table[[codon]] <- "Q"
        } else if(codon %in% c("CGT", "CGA", "CGC", "CGG", "AGA", "AGG")){
            aa_table[[codon]] <- "R"
        } else if(codon %in% c("AGC", "AGT", "TCA", "TCC", "TCG", "TCT")){
            aa_table[[codon]] <- "S"
        } else if(codon %in% c("ACA", "ACC", "ACG", "ACT")){
            aa_table[[codon]] <- "T"
        } else if(codon %in% c("GTT", "GTC", "GTA", "GTG")){
            aa_table[[codon]] <- "V"
        } else if(codon %in% c("TGG")){
            aa_table[[codon]] <- "W"
        } else if(codon %in% c("TAT", "TAC")){
            aa_table[[codon]] <- "Y"
        } else if(codon %in% c("TAA", "TAG", "TGA", "TGG")){
            aa_table[[codon]] <- "STOP"
        }
    }


    for(gene in unique(gX$sym)){

        # CALCULATE GENE-SPECIFIC BACKGROUND MUTATION RATES AND SCORES

        nuc <- c("A", "T", "C", "G")

        background_silent <- 0
        background_mis_nons <- 0

        codon_mutation_freq <- 3 * 1.1e-08 * (1 - 1.1e-08) * (1 - 1.1e-08)  # human genome mutation rate is estimated to
                                                                            # be 1.1e-08 per site per generation, according to
                                                                            # Roach et al. (2010). Science. 328 (5978): 636–9. PMID: 20220176

        # THE EXACT NUMBERS FOR NUMBER OF PEOPLE, NUMBER OF SIMULATIONS, AND LAMBDA DO NOT MATTER because
        # they are cancelled out in the missense+nonsense/silent ratio anyway. Chose these to be quite high
        # to get non-zero counts.
        num_sims = 100
        num_people = 10
        background_mut_vector <- c()

        # COMPILE GENE-SPECIFIC CODON USAGE TABLE FROM THE GENE SEQUENCE AND LENGTH INFORMATION FROM THE GENEDATA OBJECT
        # assume GeneData is structured with gene symbols as row names
        # assume GeneData has info for all of the genes in gX
        split_seq <- unlist(strsplit(GeneData[gene, "cds"], "(?<=...)", perl=TRUE))

        for(i in 1:num_sims){
            background_silent <- 0
            background_mis_nons <- 0


            for(j in 1:num_people){

                for(cod in split_seq){
                    if(rpois(1, lambda = 0.01) == 0){ # could use codon_mutation_freq as lambda (this is a closer
                                                      # reflection of real mutation rate in the human genome), but this
                                                      # results in a very low rate of incidence of mutations in the simulation.
                                                      # So, 0.01 is used instead. Again, note that the particular values are
                                                      # are not crucial, as the ratio of missense+nonsense/silent balances it out.
                        next
                    }
                    cod_nt <- strsplit(cod, "")[[1]]
                    pos <- sample(1:3, 1)
                    nucleotide <- cod_nt[pos]
                    replacement <- sample(nuc[-which(nuc == nucleotide)], 1)

                    if(pos == 1){
                        new_cod <- paste0(replacement, cod_nt[2], cod_nt[3])
                    } else if(pos == 2){
                        new_cod <- paste0(cod_nt[1], replacement, cod_nt[3])
                    } else if(pos == 3){
                        new_cod <- paste0(cod_nt[1], cod_nt[2], replacement)
                    }

                    # Assess if substitution of replacement nucleotide causes a silent mutation
                    if(aa_table[[cod]] == aa_table[[new_cod]]){
                        # If silent, add 1 to silent count
                        background_silent <- background_silent + 1
                    } else{
                        # If not silent, add 1 to missense+nonsense count
                        background_mis_nons = background_mis_nons + 1
                    }


                }

            }
            back_mut <- (background_mis_nons + 1) / (background_silent + 1) # expected rate of missense+nonsense/silent; add 1 in case of 0
            background_mut_vector <- c(background_mut_vector, back_mut)

        }

        back_mean <- mean(background_mut_vector)
        back_sd <- sd(background_mut_vector)


        # CALCULATE OBSERVED MUTATION RATE FOR GENE

        silent_count <- sum(gX$count[gX$class == "Silent" & gX$sym == gene])
        missense_count <- sum(gX$count[gX$class == "Missense_Mutation" & gX$sym == gene])
        nonsense_count <- sum(gX$count[gX$class == "Nonsense_Mutation" & gX$sym == gene])

        # observed rate of missense+nonsense/silent
        observed_mut = (missense_count + nonsense_count + 1) / (silent_count + 1) # add 1 to avoid case of division by 0



        # CALCULATE SCORE AS A LIKELIHOOD RATIO ASSUMING BACKGROUND EXPECTATION IS NORMALLY DISTRIBUTED

        # chance that this observation falls into background distribution (from simulation)
        null_prob <- pnorm(observed_mut, mean = back_mut, sd = back_sd, lower.tail=FALSE)

        # chance that this obervation falls outside the background distribution
        likelihood_ratio <- (1 - null_prob)/null_prob



        # DOES NOT MAKE SENSE TO HAVE FDR HERE, BECAUSE JUST RETURNING LIKELIHOOD RATIO. NOT DECIDING IF
        # ANYTHING IS TRUE OR NOT. FDR WOULD COME INTO PLAY WHEN DECIDING TO ACCEPT OR REJECT THE ENRICHMENT/DEPLETION.



        # Store gene and associated score in gH file

        gH <- new.env()
        if (is.null(gH[[gene]])){
            gH[[gene]] <- likelihood_ratio
        }

    }

    # Insert metadata
    attr(gH, "type")  <- hash
    saveRDS(gH, outputName)
    return(gH)


    # ==== WRITE LOG ===========================================================

    if(writeLog) {

        myTitle <- "score"

        score(gX, GeneData, outputName = "gH.rds", silent = FALSE, writeLog = TRUE)

        # Compile function call record
        myCall <- character()
        myCall[1] <- "score("
        myCall[2] <- sprintf("gX = \"%s\", ", gX)
        myCall[3] <- sprintf("GeneData = \"%s\", ", GeneData)
        myCall[4] <- sprintf("outputName = \"%s\", ", outputName)
        myCall[8] <- sprintf("silent = %s, ", as.character(silent))
        myCall[9] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")

        # indicate input object name(s)
        myInput = c("gX")
        myInput = c("GeneData")

        # Record progress information
        myNotes <- character()
        myNotes <- c(myNotes, sprintf("Built amino acid codon table."))
        myNotes <- c(myNotes, sprintf("Compile codon usage table using the GeneData resource's gene sequences."))
        myNotes <- c(myNotes, sprintf("Simulate background mutation rates per gene with 5 simulation runs in 10 individuals (per gene)."))
        myNotes <- c(myNotes, sprintf("Calculate observed mutation rates per gene."))
        myNotes <- c(myNotes, sprintf("Calculate score per gene as a likelihood ratio,
                                      assuming the background expectation is normally distributed"))
        myNotes <- c(myNotes, sprintf("Store gene score in gH hash: %s", outputName))

        # indicate output object name(s)
        myOutput = c("gH")

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
                 notes = myNotes,
                 output = myOutput
        )
    }
}
