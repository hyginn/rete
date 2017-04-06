# buildGeneDataList.R


# More details on the project task page:
# http://steipe.biochemistry.utoronto.ca/abc/students/index.php/BCB420_2017_Tasks/ImportGeneData

# Target data structure:
# list[[HGNC gene symbol]]
#    $name:        HGNC gene name
#    $UniProt:     canonical UniProt ID for this symbol
#    $Entrez:      Entrez ID
#    $assembly:    GRCh38
#    $chr:         Chromosome
#    $strand:      + or -
#    $geneStart:   Gene start coordinate
#    $geneEnd:     Gene end coordinate
#    $cdsStarts:   Vector
#    $cdsEnds:     Vector
#    $map:         Closures that return the position in coding sequence of
#                     a chromosome coordinate
#    $cds:         Coding sequence
#    $phastCons:   Conservation scores
#    $phyloP:      Conservation scores
#    $3Dmap:       data frame of 3D map information
#    $3D:          list of coordinate information

# Strategy: We identify the canonical transcript referenced for the HGNC symbol
# and get its detailed information via BioMart.

# Define request for the HGNC custom download CGI interface. Download:
#  gene symbols
#  gene names
#  locus types
#  chromosome
#  refseq ID
#  UniProt ID
#  ... only approved symbols.
HGNC_URL <- paste0("http://www.genenames.org/cgi-bin/download?",
                   "col=gd_app_sym&col=gd_app_name",
                   "&col=gd_locus_type&col=gd_pub_chrom_map",
                   "&col=gd_pub_refseq_ids&col=md_prot_id",
                   "&status=Approved&status_opt=2",
                   "&where=&order_by=gd_app_sym_sort",
                   "&format=text&limit=&submit=submit")

HGNCtable <- readr::read_delim(HGNC_URL, delim = "\t")

# Whitelist for Locus Type. For the Gene Data Table we use
# only protein data.
#
# cat(sprintf("\"%s\",\n", unique(HGNCtable$`Locus Type`)))

whitelistLT <- c("gene with protein product",
                 #"RNA, long non-coding",
                 # "pseudogene",
                 # "virus integration site",
                 # "readthrough",
                 # "phenotype only",
                 # "unknown",
                 # "region",
                 #"endogenous retrovirus",
                 # "fragile site",
                 "immunoglobulin gene",
                 # "immunoglobulin pseudogene",
                 # "transposable element",
                 # "RNA, micro",
                 # "RNA, ribosomal",
                 # "RNA, transfer",
                 # "complex locus constituent",
                 "protocadherin",
                 # "RNA, cluster",
                 # "RNA, misc",
                 # "RNA, small nuclear",
                 # "RNA, small cytoplasmic",
                 # "RNA, small nucleolar",
                 # "RNA, Y",
                 "T-cell receptor gene")
                 # "T-cell receptor pseudogene",
                 #"RNA, vault"

# nrow(HGNCtable)   # 40,928
# Drop all rows that are not in whitelist
HGNCtable <- HGNCtable[HGNCtable$`Locus Type` %in% whitelistLT, ]
# nrow(HGNCtable)   # 19,548

# x <- which(is.na(HGNCtable$`RefSeq IDs`))  # 3280
# x <- which(is.na(HGNCtable$`UniProt ID(supplied by UniProt)`) &
#                is.na(HGNCtable$`RefSeq IDs`))  # 59
# x <- which(is.na(HGNCtable$`UniProt ID(supplied by UniProt)`))  # 352

# Strategy: get data for all genes for which we have UniProt IDs
# from biomart:

if (!require(biomaRt)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("biomaRt")
    library("biomaRt")
}

myMart <- useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                # host="www.ensembl.org"
                # host = "uswest.ensembl.org"
                  host = "useast.ensembl.org"
                # host = "asia.ensembl.org"
                   )
# test connection
# getBM("ensembl_gene_id", "hgnc_symbol", "FUS", myMart)



# what filters are defined?
# (filters <- listFilters(myMart))
# (attributes <- listAttributes(myMart))

# get data - need to make query have keywords from one attribute page at
# a time
# attributes[attributes$page == "feature_page", 1:2]
# attributes[attributes$page == "sequences", 1:2]


# Function to retrieve the data for one row of the HGNC table

compileGeneDataItem <- function(HGNCrecord) {

    UniProtID <- HGNCrecord$`UniProt ID(supplied by UniProt)`

    # data from the "feature_page"
    myFeatures <- getBM(filters = "uniprotswissprot",
                        attributes = c("hgnc_symbol",
                                       "hgnc_trans_name",
                                       "uniprotswissprot",
                                       "chromosome_name",
                                       "strand",
                                       "transcript_length",
                                       "start_position",
                                       "end_position"),
                        values = UniProtID,
                        mart = myMart)

    # choose the longest transcript as a representative if there
    # is more than one transcript for this UniProt ID

    if (nrow(myFeatures) > 1) {
        iMaxLength <- which(as.numeric(myFeatures$transcript_length) ==
                                max(as.numeric(myFeatures$transcript_length),
                                    na.rm = TRUE))[1]
        myFeatures <- myFeatures[iMaxLength, ]
    }

    myHGNCtranscript <- myFeatures$hgnc_trans_name

    # data from the "sequences" page
    myCdsExons <- getBM(filters = "hgnc_trans_name",
                        attributes = c("genomic_coding_start",
                                       "genomic_coding_end"),
                        values = myHGNCtranscript,
                        mart = myMart)

    # the actual coding sequence
    myCDS<- getBM(filters = "hgnc_trans_name",
                  attributes = c("coding",
                                 "hgnc_trans_name"),
                  values = myHGNCtranscript,
                  mart = myMart)[1, 1]

    # Assemble into list

    myGDitem <- list(
        sym =         HGNCrecord$`Approved Symbol`,
        name =        HGNCrecord$`Approved Name`,
        UniProt =     UniProtID,
        RefSeq =      HGNCrecord$`RefSeq IDs`,
        assembly =   "GRCh38",
        chr =         myFeatures$chromosome_name,
        strand =      c("-", "", "+")[sign(as.numeric(myFeatures$strand)) + 2],
        geneStart =   myFeatures$start_position,
        geneEnd =     myFeatures$end_position,
        cdsStarts =   myCdsExons[ , 1],
        cdsEnds =     myCdsExons[ , 2],
        map =         NULL,
        cds =         myCDS,
        phastCons =   NULL,
        phyloP =      NULL,
        xyzMap =      NULL,
        xyz =         NULL)

    return(myGDitem)

}


# Examples:
GD1 <- compileGeneDataItem(HGNCtable[1, ])
dput(GD1, file="../A1BG.txt")
dput(compileGeneDataItem(HGNCtable[4, ]), file="../A2ML1.txt")




# [END]
