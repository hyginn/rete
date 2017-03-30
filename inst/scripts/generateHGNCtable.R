# generateHGNCtable.R

# Generates a data frame containing all valid HGNC gene symbols, calls fmatch
# on the data frame so it will create a hash of the data frame to be used by
# .fastCheckFactory() to build a validator closure for HGNC symbols.
#

# Define request for the HGNC custom download CGI interface. Download gene
# symbols and locus types only.
HGNC_URL <- paste0("http://www.genenames.org/cgi-bin/download?",
                   "col=gd_app_sym&col=gd_locus_type",
                   "&status=Approved&status_opt=2",
                   "&where=&order_by=gd_app_sym_sort",
                   "&format=text&limit=&submit=submit")

HGNCtable <- readr::read_delim(HGNC_URL, delim = "\t")

# Whitelist for Locus Type to be included. We  restrict the locus types
# to those which we consider interpreatble in the context of cancer gene
# network analysis.
#
# cat(sprintf("\"%s\",\n", unique(HGNCtable$`Locus Type`)))

whitelistLT <- c("gene with protein product",
                 "RNA, long non-coding",
                 # "pseudogene",
                 # "virus integration site",
                 # "readthrough",
                 # "phenotype only",
                 # "unknown",
                 # "region",
                 "endogenous retrovirus",
                 # "fragile site",
                 "immunoglobulin gene",
                 # "immunoglobulin pseudogene",
                 # "transposable element",
                 "RNA, micro",
                 "RNA, ribosomal",
                 "RNA, transfer",
                 # "complex locus constituent",
                 "protocadherin",
                 # "RNA, cluster",
                 # "RNA, misc",
                 "RNA, small nuclear",
                 # "RNA, small cytoplasmic",
                 "RNA, small nucleolar",
                 "RNA, Y",
                 "T-cell receptor gene",
                 # "T-cell receptor pseudogene",
                 "RNA, vault")

# Create a data frame with a single column containing the gene symbols.
HGNCsymbols <- data.frame(symbol =
        HGNCtable$`Approved Symbol`[HGNCtable$`Locus Type` %in% whitelistLT],
                          stringsAsFactors = FALSE)

# Call fmatch, so that it will attach a hash to the data frame.
fastmatch::fmatch("1", HGNCsymbols$symbol)

# save as RDS
# saveRDS(HGNCsymbols, file = "inst/extdata/HGNCsymbols.RDS")


# [END]
