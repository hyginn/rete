---
title: Rete Tutorial
bibliography: ./tutorial-draft.bib
...

# Scientific problem
- Increased understanding of which gene mutations contribute to cancer development is crucial both for genetic screening in preventive purposes and drug development to slow down disease progression.
- There are many approaches to identify single genes that are significantly mutated across different cancers [@dees_music:_2012_; @lawrence_mutational_2013; @mermel_gistic2.0_2011]. However, searching for significant mutations at the level of individual genes fails to identify functionally equivalent mutations that might occur in different genes.
    - For example, mutations to different genes could inhibit the same pathway, resulting in the pathway being altered in a high number of cancers although each individual gene is mutated at a low frequency.
    - This is especially troublesome in cancer development, where mutational heterogeneity is observed to be high across cancer samples and even between tumors of the same cancer [@vogelstein_cancer_2013; @hanahan_hallmarks_2011].
    - On the individual gene level, these rare somatic mutations which are important for cancer development, can occur with the same frequency as, and be hard to distinguish from, non-essential passenger mutations.
- To combat this, many current approaches [@printz_aacr_2017; @network_corrigendum:_2013] have tried to identify pathways of interest and look enrichment of mutations within these pathways.
    - However, cross-talk between pathways makes it difficult to assign genes strictly only to one pathway, or even one biological function.
    - Aside from incorrect pathways assignment of genes, limiting the search to specific pathways, might introduce bias which compromises the results and misses key regulatory genes. This is especially concern as current pathway annotation of genes is incomplete.

## Alternate approach
- To overcome these difficulties, our approach ([Rete](https://github.com/hyginn/rete)) is based on the methods developed by Leiserson _et al_ [@leiserson_pan-cancer_2015], which divides genes into subnetworks based on their known physical interactions with each other, rather than assigning them to any particular biological pathway. Such subnetworks can span multiple signalling pathways and render a more complete picture of what biological functions are targeted by cancerous mutations, and how such mutations can affect multiple pathways simultaneously.
- To find these commonly altered subnetworks, significantly mutated genes and their connectivity to other genes needs to be identified.
- At a high level, this process proceeds as follows
    1. Retrieve information on the mutational frequency of genes in cancer samples, and the known interactions of the protein products of these genes.
    2. Distribute the mutational frequency from a single gene across the network of genes.
    3. Identify significant subnetworks in this distributed gene network graph.

# Input data
For validation, we use the same data set used by Leiserson _et al_, consisting of ovarian adenocarcinoma samples from the TCGA [@cancer_genome_atlas_research_network_integrated_2011]. This data includes of 489 high-grade serous ovarian adenocarcinomas. Specifically, we select  one of the four subnetworks discovered by Vandin _et al._ [@vandin_discovery_2011], subnetwork S1. The genes in this network have independently been identified to overlap with genes in the focal adhesion pathway and significantly correlated with overall survival in ovarian adenocarcinoma compared to what is expected by chance [@crijns_survival-related_2009]. This data set also contains outliers, or cold nodes, as identified by Vandin _et al_.

| Reported Symbol | HGNC gene symbol | UniProt ID | Approved Name                          | Note |
|-----------------|------------------|------------|----------------------------------------|------|
| ADAM9           | ADAM9            | Q13443     | ADAM metallopeptidase domain 9         | -  |
| ITGAV           | ITGAV            | P06756     | integrin subunit alpha V               | -  |
| ITGA6           | ITGA6            | P23229     | integrin subunit alpha 6               | -  |
| ITGA3           | ITGA3            | P26006     | integrin subunit alpha 3               | -  |
| ITGB5           | ITGB5            | P18084     | integrin subunit beta 5                | -  |
| LIMK1           | LIMK1            | P53667     | LIM domain kinase 1                    | -  |
| FGFR2           | FGFR2            | P21802     | fibroblast growth factor receptor 2    | -  |
| DLST            | DLST             | P36957     | dihydrolipoamide S-succinyltransferase | -  |
| UMPS            | UMPS             | P11172     | uridine monophosphate synthetase       | -  |
| PAK4            | PAK4             | O96013     | p21 (RAC1) activated kinase 4          | -  |
| GATAD2A         | GATAD2A          | Q86YP4     | GATA zinc finger domain containing 2A  | -  |
| C2orf65         | M1AP             | Q8TC57     | meiosis 1 associated protein           | Outlier  |
| DOK1            | DOK1             | Q99704     | docking protein 1                      | Outlier  |
| DQX1            | DQX1             | Q8TE96     | DEAQ-box RNA dependent ATPase 1        | Outlier  |
| LOXL3           | LOXL3            | P58215     | lysyl oxidase like 3                   | Outlier  |
| SEMA4F          | SEMA4F           | O95754     | semaphorin 4F                          | Outlier  |

**Table 1** All the members of subnetwork 1 from Vandin _et al_.

## Mutations
- Mutations are often categorized as Single Nucleotide Variations (SNVs) or Copy Number Aberrations (CNAs). SNVs are changes changes in just a single nucleotide, while CNAs can be multiplication of anything from two base pairs to entire genes. These types of mutations are common and only a few of them are associated with cancerous outcomes.
- Data of CNAs known to be associated with different cancers can be downloaded from the [Project Genie](http://www.aacr.org/Research/Research/Pages/aacr-project-genie.aspx) and [Firehose](http://archive.broadinstitute.org/cancer/cga/Firehose) databases.
- Data of SNVs important for cancer development can be downloaded from the Project Genie, [COSMIC](http://cancer.sanger.ac.uk/cosmic) and [The Cancer Genome Database](https://cancergenome.nih.gov/).
- Examples of the input data format of the SNV and CNA data sets can be found in the [Rete GitHub repository](https://github.com/hyginn/rete/tree/master/inst/extdata).

## Protein interactions
- In addition to knowing which genes are commonly mutated in cancers, it is crucial to know how these genes can affect other genes and cellular functions.
- Most gene communication is mediated by interactions of their protein products, so this information can be approximated by our current knowledge of which proteins can bind each other and are thus could possibly interact in biologically meaningful ways.
- Protein protein interactions (PPI) can be downloaded in the form of PPI network files, for example from STRING [@szklarczyk_string_2015]and MultiNet [@khurana_interpretation_2013]].
- Examples of the input data format for the STRING PPI data can be found in the [Rete GitHub repository](https://github.com/hyginn/rete/tree/master/inst/extdata).

# Analysis
- IMPORT-M Read in mutation data of both SNVs and CNAs from their respective online repositories.
- FILTER Filter out mutations that are not associated with the cancerous phenotype, such as hypermutator genes which are known to be altered under conditions not associated with disease progression.
- COMBINE Combine the different types of mutation data into a common table containing the given mutation rate for all genes of interest.
- SCORE Assign genes a score, referred to as heat, based on their mutation frequency by scoring the combined dataset. Existing algorithms, such as MutSigCV [@lawrence_mutational_2013], can be used for this task.
- IMPORT-M Read in the PPI data from the online databases.
- ANNOTATE Annotate each vertex (protein) in the PPI network with the heat score from the corresponding gene.
- DIFFUSE Distribute the heat from one protein to its neighboring proteins according to the PPI connections. The heat distribution can be in form of a random walk to all neighbours and a probability for restarting the walk from the source. The heat flowing through an edge at each iteration is added to the total influence score of this edge. This process will highlight which paths in the PPI that mutations of a single gene can influence.
- THRESH Calculate the thresholds for FINDSUB to consider an edge for removal due to not being part of a significantly mutated path in the network.
- FINDSUB Extract the significantly mutated subnetworks from each graph by removing edges that are below the threshold calculated from THRESH. After removal of edges, subnetworks are extracted by finding the strongly connected components in the remaining network graph. Each of these networks can be assigned an aggregated heat score by summing the influence of the individual edges in the network.
- CONSENSE Find the consensus networks from all the subnetworks returned by FINDSUB, by weighing each edge by the number of subnetworks where it is present. Also identify linker nodes, which are genes that exist in multiple networks and might enable crosstalk between biological pathways.
- Downstream analyses of statistical significance and comparison to known signalling pathways by pathway enrichment scores.

# Results and interpretation
- The final results will be a list of genes constituting subnetworks and an indication of which genes facilitate crosstalk between these networks.
    - These genes have been derived from a combination the number of cancer samples that have mutations in the genes, and the interactions between genes in the subnetwork according to the PPI network.
    - As our method of edge removal uses the same heat diffusion approach as in Leiserson _et al_, we expect to arrive at similar results for the ovarian adenocarcinoma validation data set.
    - However, our implementation is flexible, in that it allows for alternate methods of edge removal to be implemented, including machine learning approaches, which could be an extension to better identify what constitute a relevant edge.
        - Machine learning approaches could be more successful in detecting edges that are relevant in driving cancer development, by integrating additional information about the connected nodes. Such node features could include the betweeness centrality, degree, semantic similarity of GO terms, and independent validation in other types of experiments, such as co-expression data sets.
        - A key challenge with machine learning approaches is to identify correctly labelled negative and positive data to use for training. For negative data, there are initiatives like the negatome which lists proteins known _not_ to interact, and positive data can be extracted from individual studies confirming the importance for disease progress of communications over a specific edge.
    - The linker genes that enable crosstalk between subnetworks might be of specific interest as key facilitators of the disease progression.
- The consensus subnetworks will contain genes that can be compared to biologically known pathways to see if they correspond closely to already charted biological functions, and if they indicate significant communication between pathways.
-  The strength of this approach is primarily to pick up networks and pathways that are commonly affected by heterogeneous cancer mutations, which could not be detected at the single gene level. Here, the overall network mutations that were previously distributed across the many genes in the network, will be aggregated and the effected of mutating any important gene in this network can be assessed.
    - Conversely, one could miss single genes that are highly important by themselves, owing to problems with the incomplete PPI networks (this gene is not well connected, and consequently cannot create significantly mutated subnetworks of large enough size to be considered in the analyses. To remedy this, once could use this approach in combination with existing single mutation scoring approaches such as MutSig and combine the results.
- There might be nodes that are not physically connected and thus not annotated in the PPI network, but which might share a functional path of edges via other nodes. Such nodes include subunits of the ribosome which not all bind to each other, but function together as a complex. If these complexes are large, the heat on each individual edge could be minor and potentially not picked through the random walk heat diffusion. Likewise, if there are many outgoing edges from the same complex, but from different subunits, these might receive a low score individually. In these cases, it could be beneficial to consider a strategy where known complexes are grouped together as one network component instead of many smaller ones. Complex data could be attained from the IntAct database [@orchard_mintact_2014].
- Where there is a hub cancer gene, like p-53, it will become extremely hot and diffuse heat to all neighbours, generating a star shaped network that is likely not of biological significance. This is largely mitigated by the restart probability of the random walk heat diffusion, but there might still be false positives with this distinct network shape, which is important to keep in mind when exploring the results.
- As always, the quality of input data is crucial for the accuracy of the suggested genes and one should be careful drawing conclusions about genes where the PPI information is known to be especially incomplete.
- Note that there are other ways of transferring information between genes than via protein to protein interactions, for example protein to gene interactions, such as is the case of transcription factor proteins binding to the DNA and regulating its expression. These interactions could be added as a future extension of the network, through the analyses of Chip-seq data.
    - Possible interaction as indicated by a PPI network also does not necessarily mean biologically meaningful interactions. Gene co-expression data, tissue/cell co-localization studies, and gene location on the genome could be used to increase the confidence that the possible interaction between two proteins would produce biologically useful information.

# References
Dees, N.D., Zhang, Q., Kandoth, C., Wendl, M.C., Schierding, W., Koboldt, D.C., Mooney, T.B., Callaway, M.B., Dooling, D., Mardis, E.R., et al. (2012). MuSiC: identifying mutational significance in cancer genomes. Genome Res. 22, 1589–1598.
Cancer Genome Atlas Research Network (2011). Integrated genomic analyses of ovarian carcinoma. Nature 474, 609–615.
Cancer Genome Atlas Research Network (2011). Integrated genomic analyses of ovarian carcinoma. Nature 474, 609–615.
Crijns, A.P.G., Fehrmann, R.S.N., de Jong, S., Gerbens, F., Meersma, G.J., Klip, H.G., Hollema, H., Hofstra, R.M.W., te Meerman, G.J., de Vries, E.G.E., et al. (2009). Survival-related profile, pathways, and transcription factors in ovarian cancer. PLoS Med. 6, e24.
Hanahan, D., and Weinberg, R.A. (2011). Hallmarks of cancer: the next generation. Cell 144, 646–674.
Khurana, E., Fu, Y., Chen, J., and Gerstein, M. (2013). Interpretation of Genomic Variants Using a Unified Biological Network Approach. PLOS Computational Biology 9, e1002886.
Lawrence, M.S., Stojanov, P., Polak, P., Kryukov, G.V., Cibulskis, K., Sivachenko, A., Carter, S.L., Stewart, C., Mermel, C.H., Roberts, S.A., et al. (2013). Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature 499, 214–218.
Leiserson, M.D.M., Vandin, F., Wu, H.-T., Dobson, J.R., Eldridge, J.V., Thomas, J.L., Papoutsaki, A., Kim, Y., Niu, B., McLellan, M., et al. (2015). Pan-cancer network analysis identifies combinations of rare somatic mutations across pathways and protein complexes. Nat Genet 47, 106–114.
Mermel, C.H., Schumacher, S.E., Hill, B., Meyerson, M.L., Beroukhim, R., and Getz, G. (2011). GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biol. 12, R41.
Network, T.C.G.A.R. (2013). Corrigendum: Comprehensive genomic characterization defines human glioblastoma genes and core pathways. Nature 494, 506.
Orchard, S., Ammari, M., Aranda, B., Breuza, L., Briganti, L., Broackes-Carter, F., Campbell, N.H., Chavali, G., Chen, C., del-Toro, N., et al. (2014). The MIntAct project--IntAct as a common curation platform for 11 molecular interaction databases. Nucleic Acids Res. 42, D358-363.
Printz, C. (2017). AACR releases large cancer genomic data set from project GENIE. Cancer 123, 1685.
Szklarczyk, D., Franceschini, A., Wyder, S., Forslund, K., Heller, D., Huerta-Cepas, J., Simonovic, M., Roth, A., Santos, A., Tsafou, K.P., et al. (2015). STRING v10: protein-protein interaction networks, integrated over the tree of life. Nucleic Acids Res. 43, D447-452.
Vandin, F., Clay, P., Upfal, E., and Raphael, B.J. (2011). Discovery of mutated subnetworks associated with clinical data in cancer. In Biocomputing 2012, (WORLD SCIENTIFIC), pp. 55–66.
Vogelstein, B., Papadopoulos, N., Velculescu, V.E., Zhou, S., Diaz, L.A., and Kinzler, K.W. (2013). Cancer genome landscapes. Science 339, 1546–1558.
