# Scientific problem
- Increased understanding of which gene mutations contribute to cancer development is crucial both for genetic screening in preventive purposes and drug development to slow down disease progression.
- There are many approaches to identify single genes that are significantly mutated across different cancers [examples, refs]. However, searching for significant mutations at the level of individual genes misses fails to identify functionally equivalent mutations that might occur in different genes.
    - For example, mutations to different genes could inhibit the same pathway, resulting in the pathway being altered in a high number of cancers although each individual gene is mutated at a low frequency.
    - This is especially troublesome in cancer development, where mutational heterogeneity is observed to be high across cancer samples and even between tumors of the same cancer [ref].
    - On the individual gene level, these rare somatic mutations which are important for cancer development, can occur with the same frequency as, and be hard to distinguish from, non-essential passenger mutations.
- To combat this, many current approaches [ref] have tried to identify pathways of interest and look enrichment of mutations within these pathways.
    - However, cross-talk between pathways makes it difficult to assign genes strictly only to one pathway, or even one biological function.
    - Aside from incorrect pathways assignment of genes, limiting the search to specific pathways, might introduce bias which compromises the results and misses key regulatory genes. This is especially concern as current pathway annotation of genes is incomplete.

## Our approach
- To overcome these difficulties, our approach divides genes into subnetworks based on their known physical interactions with other genes, rather than assigning them to any particular biological pathway. Such subnetworks can span multiple pathways and give a clearer picture of what biological functions are targeted by cancerous mutations.
- To find these commonly mutated subnetworks, significantly mutated genes and their connectivity to other genes needs to be identified.
- At a high level this process looks like this.
    1. Retrieve information on the mutational frequency of genes in cancer samples, and the known interactions of the protein products of these genes.
    2. Distribute the mutational frequency from a single gene across the network of genes.
    3. Identify significant subnetworks in this distributed gene network graph.

# Input data

## Mutations
- Mutations are often categorized as Single Nucleotide Variations (SNVs) or Copy Number Aberrations (CNAs). SNVs are changes changes in just a single nucleotide, while CNAs can be multiplication of anything from two base pairs to entire genes. These types of mutations are common and only a few of them are associated with cancerous outcomes.
- Data of CNAs known to be associated with different cancers can be downloaded from the Project Genie [ref] and Firehose [ref] databases.
- Data of SNVs important for cancer development can be downloaded from the Project Genie, COSMIC and The Cancer Genome Database [refs].

## Protein interactions
- In addition to knowing which genes are commonly mutated in cancers, it is crucial to know how these genes can affect other genes and cellular functions.
- Most gene communication is mediated by interactions of their protein products, so this information can be approximated by our current knowledge of which proteins can bind each other and are thus could possibly interact in biologically meaningful ways.
- Protein protein interactions (PPI) can be downloaded in the form of PPI network files, for example from STRING, MultiNet, and MITAB [refs].

# Analysis
- IMPORT-M Read in mutation data of both SNVs and CNAs.
- FILTER Filter out mutations that are not associated with the cancerous phenotype, such as hypermutator genes.
- COMBINE Combine the different types of mutation data into one table containing the given mutation rate for all genes of interest.
- SCORE Assign genes a score, referred to as heat, based on their mutation frequency by scoring the combined dataset. Can use existing algorithms, such as MutSigCV [ref].
- IMPORT-M Read in the PPI data.
- ANNOTATE each vertex (protein) in the PPI network with the heat scores from the corresponding gene.
- DIFFUSE Distribute the heat from one protein to its neighboring proteins according to the PPI connections. The heat distribution can be in form of a random walk to all neighbours and a probability for restarting the walk from the source. This will highlight which paths in the PPI that mutations to a single gene can influence.
- THRESH Calculate the thresholds to use for FINDSUB to consider an edge for removal due to not being part of a significantly mutated path in the network.
- FINDSUB Extract the significantly mutated subnetworks from each graph, by removing edges that are below the threshold calculated from THRESH.
- CONSENSE Find the consensus subnetwork from all the subnetworks, by weighing each edge by the number of subnetworks where it is present.
- Downstream analyses of statistical significance and comparison to known pathways by pathway enrichment scores.

# Results and interpretation
- The final results will be a list of genes constituting subnetworks and an indication of which genes facilitate crosstalk between these networks.
    - The linker genes that enable crosstalk between subnetworks might be of specific interest as key facilitators of the disease progression.
- These genes have been derived from a combination the number of cancer samples that have mutations in the genes, and the interactions between genes in the subnetwork according to the PPI network.
- These genes can be compared to biologically known pathways to see if they correspond closely to already charted biological function, and if they indicate significant communication between pathways.
- Where there is a hub cancer gene, like p-53, it will become extremely hot and diffuse heat to all neighbours, generating a star shaped network that is likely not of biological significance. This is largely mitigated by the restart probability of the random walk heat diffusion, but there might still be false positives with this distinct network shape, which is important to keep in mind when exploring the results.
-  The strength of this approach is primarily to pick up networks and pathways that are commonly affected by heterogeneous cancer mutations, which could not be detected at the single gene level. Here, the overall network mutations that were previously distributed across the many genes in the network, will be aggregated and the effected of mutating any important gene in this network can be assessed.
    - Conversely, one could miss single genes that are highly important by themselves, owing to problems with the incomplete PPI networks (this gene is not well connected, and consequently cannot create significantly mutated subnetworks of large enough size to be considered in the analyses. To remedy this, once could use this approach in combination with existing single mutation scoring approaches such as MutSig and combine the results.
- As always, the quality of input data is crucial for the accuracy of the suggested genes and one should be careful drawing conclusions about genes where the PPI information is known to be especially incomplete.
-  There might be nodes that are not physically connected and thus not annotated in the PPI network, but which might share a functional path of edges via other nodes. Such nodes include subunits of the ribosome which not all bind to each other, but function together as a complex. If these complexes are large, the heat on each individual edge could be minor and potentially not picked through the random walk heat diffusion. Likewise, if there are many outgoing edges from the same complex, but from different subunits, these might receive a low score individually. In these cases, it could be beneficial to consider a strategy where known complexes are grouped together as one network component instead of many smaller ones.
- Note that there are other ways of transferring information between genes than via protein to protein interactions, for example protein to gene interactions, such as is the case of transcription factor proteins binding to the DNA and regulating its expression. These interactions could be added as a future extension of the network, through the analyses of Chip-seq data.
    - Possible interaction as indicated by a PPI network also does not necessarily mean biologically meaningful interactions. Gene co-expression data, tissue/cell co-localization studies, and gene location on the genome could be used to increase the confidence that the possible interaction between two proteins would produce biologically useful information.

