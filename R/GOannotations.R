#' GOannotation.R
#'
#' Produced annotated GO terms graph with provided GO terms and GO annotations
#'
#' \code{GOannotation} Builds a directed acyclic graph from inputs of data files for
#' GO terms and GO annotations and save two data frame files. The first dataframe is the DAG built
#' based on input and the second dataframe is genes, which appeared in the DAG, with associated
#' GO terms that they appear in
#'
#' @param fNameGO A vector of path of fully qualifed filenames of GO terms .obo files.
#' @param fNameGOA A vector of path of fully qualifed filenames of GO annotation .gaf files.
#' @param fGDAG The fully qualified filename/path of a GO DAG output file, by default is GDAG.RDS.
#' @param fGg The fully qualified filename/path of a genes with associated GO terms output file, by default is Gg.RDS.
#' @param silent Boolean option for writing graph building process information to console, FALSE by default.
#' @param writeLog Boolean option for log results, TRUE by default.
#'
#' @return The size of each data frame. The number of nodes and storage requirement for the first dataframe.
#' The number of genes and storage requirements for the second dataframe.
#'
#' @examples
#' \dontrun{GOannotations(c("go.obo"), c("goa_human_complex.gaf")}

#' @export
GOannotation <- function(fNameGO, fNameGOA, FGDAG="GDAG.RDS", FGg="Gg.RDS", silent=FALSE, writeLog=TRUE) {
    # Download Ontology (GO) data
    # download.file("http://purl.obolibrary.org/obo/go.obo", destfile = "inst/extdata/go.obo")

    # Download GOA data (GOA)
    # ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/

    # Pseudocode begins
    # Packages for using readGAF, GOA files
    # source("https://bioconductor.org/biocLite.R")
    # biocLite("mgsa")
    # biocLite("GO.db")

    # GOA <- readGAF(fNameGOA) store all annotations in a object
    # create a hashtable, hashTable <- hash::hash()

    # Parse obo files
    # obo <- file("./inst/extdata/go.obo", open="r")
    # lines <- readlines(obo)
    # number of GO terms is sum(l=="[Term]")
    # getting indices of each term in lines, terms <- which(l=="[Term]")
    # iterate through terms, by increasing index in lines:
    #     get GO id from keyword "id" and use it as hash for hashtable
    #
    #     make a vector with name parent from keyword "is_a"
    #
    #     get gene symbols from GOA
    #     sample code:
    #     rawSymbol <- paste(GOA@itemAnnotations[GOA@sets[names(GOA@sets)==<current GO id>][[1]],]$symbol)
    #     for (gene in rawSymbol) {
    #         geneSymbol <- strsplit(gene, "_")[[1]][1]
    #         if isHGNCsymbole(geneSymbol) is TRUE:
    #             append the geneSymbol to a vector with name gene symbol
    #         otherwise, add geneSymbol to a vector called nonHGNC
    #     }
    #
    #     get the current GO term name using keyword "name"
    #
    #
    #     if blank line encountered:
    #         assign the above three information to the hash as a list
    #         jump to next idex for "[Term]"

    # infer children list
    # for each key in the hashtable:
    #     add its own GO id with vector named "name" to the parent's children vector
    #     check if the terms with empty parent vector are three, which are GO:0003674, GO:0005575, and GO:0008150
    #     report error if the check failed

    # compute distance and prune leaves
    # initialize a vector for leaves
    # using the helper function below to enqueue
    # for each root:
    #     enqueue(root)
    #     while the queue is not empty:
    #       if the length of current node's gene symbol vector is 0:
    #           if the node's children and parents vector are exactly 1:
    #               add its name to its' children name in the form <GO:xxxxx:GO:yyyyy> and
    #               change its name accordingly in its parent's children vector.
    #       if current node's parent vector is empty: set the minumum distance to 0
    #       otherwise, set the minimum distance to min(all parents' minimum distance) + 1
    #
    #       if the current node's children vector is empty, append its key to the vector for leaves
    #       dequeue(current node)
    #       enqueue(current node$children)

    # compute gene counts
    # get a vector of hash keys which are leave from previous BFS
    # for every gene symbol in the leaves vector:
    #       set the current node's gene count to the sum of all its children gene count + length of its gene vector
    #       if the symbol has not appeared in a vector called unique_genes
    #           add the symbol to the vector
    #       dequeue the current node
    #       delete current node from its parents' children vector
    #       if any parent's children vector length is one:
    #           enqueue the parent

    # building the second hash
    # initialize another hashtable2
    # for every gene in the unique gene vector obtained above:
    #       use the gene symbol as hash key and push into hashtable2
    #       create a list of 3 vectors with names GO:0003674, GO:0005575, and GO:0008150 as attributes to hash key
    #       for each root:
    #           go down the DAG using top down BFS to check whether the gene symbol is in the gene symbol vector
    #           record the GO term id to the corresponding vector

    # create two data frames from the two hashtables
    # for DAG, the rows will be GO id terms and columns will be name, gene symbols,
    #           parents, children, gene numbers and minimum distance to root
    # for the gene hashtable, the rows will be genes and the columns will be GO id terms, so the value will be 1
    #           if the gene is annotated to such term, -1 otherwise.

    # return the length of GO term ids and length of the genes, in other words, the length of each hashtable
    # as a tuple, e.g. (500, 300)

    # write to log file
    # save RDS files to path of FGDAG and FGg

}

# helper_function for queue structure
# inspired by
# enqueue <- function(queue, new.element) {
#     temp.queue <- new.env()
#     temp.queue$value <- new.element
#     temp.queue$next.element <- queue$next.element
#     queue$next.element <- temp.queue
#     queue$value <- NULL
# }
#
# dequeue <- function(queue) {
#     value <- queue$next.element$value
#     queue$next.element <- (queue$next.element)$next.element
#     queue$value <- NULL
#     return(value)
# }

# [END]
