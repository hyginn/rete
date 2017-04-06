.compileCUT <- function(geneData) {

    # Generate codons
    aminoHash <- readRDS(system.file("extdata",
                        "aminoAcids.rds",
                        package="rete"))
    RNA <- c("A", "C", "G", "T")
    codons <- character(length = 64)
    curr <- 1
    for (first in RNA) {
        for (second in RNA) {
            for (third in RNA) {
                combo <- paste0(first, second, third)
                codons[curr] <- combo
                curr <- curr + 1
            }
        }
    }

    # Generate the CUT table
    CUT <- list()
    for (codon in codons) {
        CUT[[codon]] <- structure(list(
            A = aminoHash[[codon]],
            counts = 0,
            weights = 0.0
        ))
    }

    for (entry in geneData) {
        #if (!.isHGNCSymbol(entry$sym)) {
        #    next
        #} else {
        cds <- entry$cds
        cdsSplit <- unlist(strsplit(cds, "(?<=...)", perl=TRUE))
        result <- table(cdsSplit)
        for (codon in rownames(result)) {
            count <- result[[codon]]
            CUT[[codon]]$count <- CUT[[codon]]$count + count
        }
    }


    return(CUT)
}
