.compileCUT <- function() {
    # Generate something to associate a codon to an amino acid
    # Was thinking of manually generating a hash table with keys of codons and
    # Values of the amino acid.

    # Should also have a vector of codons

    # Generate codon structure
    CUT <- list()
    for (codon in codonList) {
        CUT[[codon]] <- structure(list(
            A = lookUp from hashTable
            # initially zero
            counts = 0
            weights = 0.0
        ))
    }

    # Metadata setup
    return(CUT)
}

.populateCUT <- function(cutObject, geneData) {
    for (gene in geneData) {
        if (!.isHGNC(gene$Name)) {
            break
        }
        cds <- gene$cds
        split <- Into Threes
        results <- table(split)
        for (codon in rownames(result)) {
            count <- results[[codon]]
            cutObject[[codon]]$count <- cutObject[[codon]]$count + count
        }

    }

}
